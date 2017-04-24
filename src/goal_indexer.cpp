#include <apf.h>
#include <apfAlbany.h>
#include <apfMesh2.h>
#include <apfMixedNumbering.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"

namespace goal {

using Teuchos::rcp;

static void destroy_local_numberings(std::vector<apf::Numbering*>& n) {
  for (size_t i = 0; i < n.size(); ++i)
    if (n[i]) apf::destroyNumbering(n[i]);
  n.resize(0);
}

static void destroy_global_numberings(std::vector<apf::GlobalNumbering*>& n) {
  for (size_t i = 0; i < n.size(); ++i)
    if (n[i]) apf::destroyGlobalNumbering(n[i]);
  n.resize(0);
}

static void synchronize_global_numberings(
    std::vector<apf::GlobalNumbering*>& n) {
  for (size_t i = 0; i < n.size(); ++i) apf::synchronize(n[i]);
}

Indexer::Indexer(RCP<Discretization> d, std::vector<RCP<Field> > f) {
  auto t0 = time();
  disc = d;
  fields = f;
  assert(fields.size() > 0);
  for (size_t i = 0; i < fields.size(); ++i) {
    assert(fields[i]->get_associated_dof_idx() == int(i));
    apf_fields.push_back(fields[i]->get_apf_field());
  }
  comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  compute_owned_map();
  compute_ghost_map();
  compute_graphs();
  compute_node_map();
  compute_coords();
  compute_node_sets();
  set_elem_block(0);
  auto t1 = time();
  print(" > indexer built in %f seconds", t1 - t0);
}

Indexer::~Indexer() {
  destroy_local_numberings(owned_numberings);
  destroy_local_numberings(ghost_numberings);
  destroy_global_numberings(global_numberings);
  owned_map = Teuchos::null;
  ghost_map = Teuchos::null;
  owned_graph = Teuchos::null;
  ghost_graph = Teuchos::null;
  node_map = Teuchos::null;
  coords = Teuchos::null;
}

void Indexer::set_elem_block(const int idx) {
  elem_block_idx = idx;
  for (int i = 0; i < get_num_dof_fields(); ++i) fields[i]->set_elem_block(idx);
}

int Indexer::get_num_dof_fields() const { return fields.size(); }

int Indexer::get_num_total_elem_dofs() const {
  int dofs = 0;
  for (int f = 0; f < get_num_dof_fields(); ++f)
    dofs += fields[f]->get_num_elem_dofs();
  return dofs;
}

int Indexer::get_dof_idx(std::string const& name) const {
  int idx = -1;
  for (int f = 0; f < get_num_dof_fields(); ++f)
    if (apf::getName(apf_fields[f]) == name) idx = f;
  if (idx < 0) fail("could not find dof field %s", name.c_str());
  return idx;
}

int Indexer::get_elem_dof_offset(const int idx, const int n, const int c) {
  int offset = 0;
  int nc = fields[idx]->get_num_components();
  for (int f = 0; f < idx; ++f) offset += fields[f]->get_num_elem_dofs();
  offset += n * nc + c;
  return offset;
}

LO Indexer::get_ghost_lid(
    const int idx, apf::MeshEntity* e, const int n, const int c) {
  static apf::NewArray<int> lids;
  int nc = fields[idx]->get_num_components();
  auto numbering = ghost_numberings[idx];
  apf::getElementNumbers(numbering, e, lids);
  return lids[n * nc + c];
}

LO Indexer::get_owned_lid(const int idx, apf::Node const& n, const int c) {
  auto numbering = owned_numberings[idx];
  return apf::getNumber(numbering, n.entity, n.node, c);
}

void Indexer::get_ghost_lids(apf::MeshEntity* e, std::vector<LO>& lids) {
  apf::getElementNumbers(ghost_numberings, e, lids);
}

int Indexer::get_num_node_sets() const { return disc->get_num_node_sets(); }

std::vector<apf::Node> const& Indexer::get_node_set_nodes(
    std::string const& set, const int idx) {
  if (!node_sets.count(set)) fail("node set %s not found", set.c_str());
  return node_sets[set][idx];
}

apf::Mesh* Indexer::get_apf_mesh() { return disc->get_apf_mesh(); }

apf::Numbering* Indexer::get_apf_numbering(const int idx) {
  assert(idx < get_num_dof_fields());
  return owned_numberings[idx];
}

apf::Numbering* Indexer::get_ghost_apf_numbering(const int idx) {
  assert(idx < get_num_dof_fields());
  return ghost_numberings[idx];
}

void Indexer::compute_owned_map() {
  destroy_local_numberings(owned_numberings);
  destroy_global_numberings(global_numberings);
  int num_owned = apf::numberOwned(apf_fields, owned_numberings);
  apf::makeGlobal(owned_numberings, global_numberings);
  Teuchos::Array<GO> indices(num_owned);
  apf::DynamicArray<apf::Node> owned_nodes;
  apf::DynamicArray<apf::Node> global_nodes;
  for (int f = 0; f < get_num_dof_fields(); ++f) {
    auto owned_numbering = owned_numberings[f];
    auto global_numbering = global_numberings[f];
    apf::getNodes(owned_numbering, owned_nodes);
    apf::getNodes(global_numbering, global_nodes);
    assert(owned_nodes.size() == global_nodes.size());
    int nc = apf::countComponents(owned_numbering);
    for (size_t n = 0; n < owned_nodes.size(); ++n) {
      auto n_lo = owned_nodes[n];
      auto n_go = global_nodes[n];
      for (int c = 0; c < nc; ++c) {
        auto ent = n_lo.entity;
        auto node = n_lo.node;
        LO owned = apf::getNumber(owned_numbering, ent, node, c);
        GO global = apf::getNumber(global_numbering, n_go, c);
        indices[owned] = global;
      }
    }
  }
  owned_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
  synchronize_global_numberings(global_numberings);
}

void Indexer::compute_ghost_map() {
  destroy_local_numberings(ghost_numberings);
  int num_ghost = apf::numberGhost(apf_fields, ghost_numberings);
  Teuchos::Array<GO> indices(num_ghost);
  apf::DynamicArray<apf::Node> ghost_nodes;
  apf::DynamicArray<apf::Node> global_nodes;
  for (int f = 0; f < get_num_dof_fields(); ++f) {
    auto ghost_numbering = ghost_numberings[f];
    auto global_numbering = global_numberings[f];
    apf::getNodes(ghost_numbering, ghost_nodes);
    apf::getNodes(global_numbering, global_nodes);
    assert(ghost_nodes.size() == global_nodes.size());
    int nc = apf::countComponents(ghost_numbering);
    for (size_t n = 0; n < ghost_nodes.size(); ++n) {
      auto n_lo = ghost_nodes[n];
      auto n_go = global_nodes[n];
      for (int c = 0; c < nc; ++c) {
        auto ent = n_lo.entity;
        auto node = n_lo.node;
        LO ghost = apf::getNumber(ghost_numbering, ent, node, c);
        GO global = apf::getNumber(global_numbering, n_go, c);
        indices[ghost] = global;
      }
    }
  }
  ghost_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
}

void Indexer::compute_graphs() {
  auto m = disc->get_apf_mesh();
  auto estimate = get_num_total_elem_dofs();
  owned_graph = rcp(new Graph(owned_map, estimate));
  ghost_graph = rcp(new Graph(ghost_map, estimate));
  std::vector<long> numbers;
  std::vector<GO> go_numbers;
  apf::MeshEntity* elem = 0;
  auto elems = m->begin(disc->get_num_dims());
  while ((elem = m->iterate(elems))) {
    apf::getElementNumbers(global_numberings, elem, numbers);
    auto ndofs = numbers.size();
    go_numbers.resize(ndofs);
    for (size_t dof = 0; dof < ndofs; ++dof) go_numbers[dof] = numbers[dof];
    for (size_t dof = 0; dof < ndofs; ++dof) {
      GO row = numbers[dof];
      auto cols = Teuchos::arrayView<GO>(&go_numbers[0], ndofs);
      ghost_graph->insertGlobalIndices(row, cols);
    }
  }
  m->end(elems);
  ghost_graph()->fillComplete();
  auto exporter = rcp(new Export(ghost_map, owned_map));
  owned_graph->doExport(*ghost_graph, *exporter, Tpetra::INSERT);
  owned_graph->fillComplete();
}

static bool should_skip(std::vector<RCP<Field> > const& f) {
  bool skip = false;
  for (size_t i = 0; i < f.size()-1; ++i) {
    if (f[i]->get_basis_type() != f[i+1]->get_basis_type()) skip = true;
    if (f[i]->get_p_order() != f[i]->get_p_order()) skip = true;
  }
  return skip;
}

static int get_num_dofs(std::vector<RCP<Field> > const& f) {
  int ndofs = 0;
  for (size_t i = 0; i < f.size(); ++i)
    ndofs += f[i]->get_num_components();
  return ndofs;
}

static int get_num_owned_nodes(std::vector<apf::Numbering*> const& n) {
  int owned = 0;
  for (size_t i = 0; i < n.size(); ++i)
    owned += apf::countNodes(n[i]);
  return owned;
}

static int get_offset(std::vector<RCP<Field> > const& f, int idx) {
  int offset = 0;
  for (int i = 0; i < idx; ++i)
    offset += f[i]->get_num_components();
  return offset;
}

static LO map_dof_to_node(
    apf::Numbering* on, apf::Node const& n, int c, int ndofs) {
  return (apf::getNumber(on, n.entity, n.node, 0)-c)/ndofs;
}

static GO map_dof_to_node(
    apf::GlobalNumbering* gn, apf::Node const& n, int c, int ndofs) {
  return (apf::getNumber(gn, n.entity, n.node, 0)-c)/ndofs;
}

void Indexer::compute_node_map() {
  if (should_skip(fields)) {
    destroy_global_numberings(global_numberings);
    return;
  }
  auto nf = get_num_dof_fields();
  auto ndofs = get_num_dofs(fields);
  auto non = get_num_owned_nodes(owned_numberings);
  Teuchos::Array<GO> indices(non);
  apf::DynamicArray<apf::Node> owned;
  for (int f = 0; f < nf; ++f) {
    auto on = owned_numberings[f];
    auto gn = global_numberings[f];
    auto c = get_offset(fields, f);
    apf::getNodes(on, owned);
    for (size_t n = 0; n < owned.size(); ++n) {
      LO nlid = map_dof_to_node(on, owned[n], c, ndofs);
      GO ngid = map_dof_to_node(gn, owned[n], c, ndofs);
      indices[nlid] = ngid;
    }
  }
  node_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
  destroy_global_numberings(global_numberings);
}

static void get_coord(apf::Field* f, apf::Node const& n, apf::Vector3& x) {
  apf::Vector3 xi(0, 0, 0);
  auto m = apf::getMesh(f);
  auto me = apf::createMeshElement(m, n.entity);
  auto e = apf::createElement(f, me);
  auto type = m->getType(n.entity);
  apf::getShape(f)->getNodeXi(type, n.node, xi);
  apf::mapLocalToGlobal(me, xi, x);
  apf::destroyElement(e);
  apf::destroyMeshElement(me);
}

void Indexer::compute_coords() {
  if (node_map == Teuchos::null) return;
  auto dim = disc->get_num_dims();
  auto nf = get_num_dof_fields();
  auto ndofs = get_num_dofs(fields);
  coords = rcp(new MultiVector(node_map, dim, false));
  apf::Vector3 point(0, 0, 0);
  apf::DynamicArray<apf::Node> owned;
  for (int f = 0; f < nf; ++f) {
    auto on = owned_numberings[f];
    auto c = get_offset(fields, f);
    apf::getNodes(on, owned);
    for (size_t n = 0; n < owned.size(); ++n) {
      auto node = owned[n];
      get_coord(apf_fields[f], node, point);
      LO nlid = map_dof_to_node(on, node, c, ndofs);
      for (int d = 0; d < dim; ++d)
        coords->replaceLocalValue(nlid, d, point[d]);
    }
  }

  MM_Writer::writeDenseFile("COORDS", coords);
}

void Indexer::compute_node_sets() {
  auto m = disc->get_apf_mesh();
  auto sets = disc->get_model_sets();
  auto nns = sets->models[0].size();
  auto nf = get_num_dof_fields();
  for (size_t s = 0; s < nns; ++s) {
    auto name = disc->get_node_set_name(s);
    node_sets[name].resize(nf);
    for (int f = 0; f < nf; ++f) node_sets[name][f].resize(0);
  }
  apf::DynamicArray<apf::Node> owned;
  for (int f = 0; f < nf; ++f) {
    auto on = owned_numberings[f];
    apf::getNodes(on, owned);
    for (size_t n = 0; n < owned.size(); ++n) {
      auto node = owned[n];
      auto ent = node.entity;
      std::set<apf::StkModel*> mset;
      apf::collectEntityModels(m, sets->invMaps[0], m->toModel(ent), mset);
      if (mset.empty()) continue;
      APF_ITERATE(std::set<apf::StkModel*>, mset, mit) {
        auto ns = *mit;
        auto nsn = ns->stkName;
        node_sets[nsn][f].push_back(node);
      }
    }
  }
}

}  /* namespace goal */
