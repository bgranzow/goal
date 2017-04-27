#include <apfAlbany.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_strided_indexer.hpp"

namespace goal {

using Teuchos::rcp;

static void validate_fields(std::vector<RCP<Field> >& f) {
  for (size_t i = 0; i < f.size()-1; ++i) {
    assert(f[i]->get_basis_type() == f[i+1]->get_basis_type());
    assert(f[i]->get_p_order() == f[i+1]->get_p_order());
  }
}

static int get_num_eqs(std::vector<RCP<Field> > const& f) {
  int num_eqs = 0;
  for (size_t i = 0; i < f.size(); ++i)
    num_eqs += f[i]->get_num_components();
  return num_eqs;
}

static void destroy_local(apf::Numbering** n) {
  if (! *n) return;
  apf::destroyNumbering(*n);
  *n = 0;
}

static void destroy_global(apf::GlobalNumbering** n) {
  if (! *n) return;
  apf::destroyGlobalNumbering(*n);
  *n = 0;
}

StridedIndexer::StridedIndexer(
    RCP<Discretization> d, std::vector<RCP<Field> > f)
    : Indexer(d, f),
      owned_numbering(0),
      ghost_numbering(0),
      global_numbering(0) {
  auto t0 = time();
  validate_fields(fields);
  num_eqs = get_num_eqs(fields);
  mesh = disc->get_apf_mesh();
  shape = apf::getShape(fields[0]->get_apf_field());
  compute_owned_maps();
  compute_ghost_map();
  compute_graphs();
  compute_coords();
  compute_node_sets();
  auto t1 = time();
  print(" > strided indexer: number pde equations: %d", num_eqs);
  print(" > strided indexer: built in %f seconds", t1 - t0);
}

StridedIndexer::~StridedIndexer() {
  destroy_local(&owned_numbering);
  destroy_local(&ghost_numbering);
  destroy_global(&global_numbering);
}

static LO get_dof(const LO node_id, const int eq, const int num_eqs) {
  return node_id*num_eqs + eq;
}

static GO get_gdof(const GO node_id, const int eq, const int num_eqs) {
  return node_id*num_eqs + eq;
}

static int get_eq_offset(
   std::vector<RCP<Field> > const& f, const int idx, const int c) {
  int eq = 0;
  for (int i = 0; i < idx; ++i)
    eq += f[i]->get_num_components();
  return eq + c;
}

int StridedIndexer::get_elem_dof_offset(
    const int idx, const int n, const int c) {
  int eq = get_eq_offset(fields, idx, c);
  return n*num_eqs + eq;
}

LO StridedIndexer::get_ghost_lid(
    const int idx, apf::MeshEntity* e, const int n, const int c) {
  int eq = get_eq_offset(fields, idx, c);
  static apf::NewArray<int> node_ids;
  apf::getElementNumbers(ghost_numbering, e, node_ids);
  return get_dof(node_ids[n], eq, num_eqs);
}

LO StridedIndexer::get_owned_lid(
    const int idx, apf::Node const& n, const int c) {
  int eq = get_eq_offset(fields, idx, c);
  LO node_id = apf::getNumber(owned_numbering, n.entity, n.node, 0);
  return get_dof(node_id, eq, num_eqs);
}

void StridedIndexer::get_ghost_lids(
    apf::MeshEntity* e, std::vector<LO>& lids) {
  lids.resize(get_num_total_elem_dofs());
  static apf::NewArray<int> node_ids;
  apf::getElementNumbers(ghost_numbering, e, node_ids);
  int elem_dof = 0;
  for (int n = 0; n < fields[0]->get_num_elem_nodes(); ++n) {
    for (int eq = 0; eq < num_eqs; ++eq) {
      lids[elem_dof] = get_dof(node_ids[n], eq, num_eqs);
      ++elem_dof;
    }
  }
}

void StridedIndexer::add_to_fields(
    std::vector<RCP<Field> > const& flds, RCP<Vector> du) {
  assert((int)flds.size() == get_num_dof_fields());
  auto data = du->get1dView();
  double values[3] = {0.0, 0.0, 0.0};
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(owned_numbering, owned);
  for (size_t n = 0; n < owned.size(); ++n) {
    auto node = owned[n];
    auto ent = node.entity;
    auto lnode = node.node;
    for (size_t f = 0; f < flds.size(); ++f) {
      auto fld = flds[f];
      auto apf_fld = fld->get_apf_field();
      auto nc = fld->get_num_components();
      apf::getComponents(apf_fld, ent, lnode, values);
      for (int c = 0; c < nc; ++c) {
        LO row = get_owned_lid(f, node, c);
        values[c] += data[row];
      }
      apf::setComponents(apf_fld, ent, lnode, values);
    }
  }
  for (size_t f = 0; f < flds.size(); ++f)
    apf::synchronize(flds[f]->get_apf_field());
}

void StridedIndexer::set_to_fields(
    std::vector<RCP<Field> > const& flds, RCP<Vector> x) {
  assert((int)flds.size() == get_num_dof_fields());
  auto data = x->get1dView();
  double values[3] = {0.0, 0.0, 0.0};
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(owned_numbering, owned);
  for (size_t n = 0; n < owned.size(); ++n) {
    auto node = owned[n];
    auto ent = node.entity;
    auto lnode = node.node;
    for (size_t f = 0; f < flds.size(); ++f) {
      auto fld = flds[f];
      auto apf_fld = fld->get_apf_field();
      auto nc = fld->get_num_components();
      for (int c = 0; c < nc; ++c) {
        LO row = get_owned_lid(f, node, c);
        values[c] = data[row];
      }
      apf::setComponents(apf_fld, ent, lnode, values);
    }
  }
  for (size_t f = 0; f < flds.size(); ++f)
    apf::synchronize(flds[f]->get_apf_field());
}

void StridedIndexer::compute_owned_maps() {
  destroy_local(&owned_numbering);
  destroy_global(&global_numbering);
  owned_numbering = apf::numberOwnedNodes(mesh, "owned", shape);
  global_numbering = apf::makeGlobal(owned_numbering, false);
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(global_numbering, owned);
  auto num_owned = owned.getSize();
  Teuchos::Array<GO> indices(num_owned);
  for (size_t n = 0; n < num_owned; ++n)
    indices[n] = apf::getNumber(global_numbering, owned[n]);
  node_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
  indices.resize(num_eqs * num_owned);
  for (size_t n = 0; n < num_owned; ++n) {
    GO gid = apf::getNumber(global_numbering, owned[n]);
    for (int eq = 0; eq < num_eqs; ++eq)
      indices[get_dof(n, eq, num_eqs)] = get_gdof(gid, eq, num_eqs);
  }
  owned_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
  apf::synchronize(global_numbering);
}

void StridedIndexer::compute_ghost_map() {
  destroy_local(&ghost_numbering);
  ghost_numbering = apf::numberOverlapNodes(mesh, "ghost", shape);
  apf::DynamicArray<apf::Node> ghost;
  apf::getNodes(ghost_numbering, ghost);
  auto num_ghost = ghost.getSize();
  Teuchos::Array<GO> indices(num_eqs * num_ghost);
  for (size_t n = 0; n < num_ghost; ++n) {
    GO gid = apf::getNumber(global_numbering, ghost[n]);
    for (int eq = 0; eq < num_eqs; ++eq)
      indices[get_dof(n, eq, num_eqs)] = get_gdof(gid, eq, num_eqs);
  }
  ghost_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
}

void StridedIndexer::compute_graphs() {
  auto dim = disc->get_num_dims();
  auto estimate = get_num_total_elem_dofs();
  owned_graph = rcp(new Graph(owned_map, estimate));
  ghost_graph = rcp(new Graph(ghost_map, estimate));
  apf::MeshEntity* elem;
  auto elems = mesh->begin(dim);
  while ((elem = mesh->iterate(elems))) {
    apf::NewArray<long> gids;
    int nnodes = apf::getElementNumbers(global_numbering, elem, gids);
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < num_eqs; ++j) {
        GO row = get_gdof(gids[i], j, num_eqs);
        for (int k = 0; k < nnodes; ++k) {
          for (int l = 0; l < num_eqs; ++l) {
            GO col = get_gdof(gids[k], l, num_eqs);
            auto col_av = Teuchos::arrayView(&col, 1);
            ghost_graph->insertGlobalIndices(row, col_av);
  }}}}}
  mesh->end(elems);
  ghost_graph->fillComplete();
  auto exporter = rcp(new Export(ghost_map, owned_map));
  owned_graph->doExport(*ghost_graph, *exporter, Tpetra::INSERT);
  owned_graph->fillComplete();
  destroy_global(&global_numbering);
}

static void get_coord(apf::Field* f, apf::Node const& n, apf::Vector3& x) {
  apf::Vector3 xi(0,0,0);
  auto m = apf::getMesh(f);
  auto me = apf::createMeshElement(m, n.entity);
  auto e = apf::createElement(f, me);
  auto type = m->getType(n.entity);
  apf::getShape(f)->getNodeXi(type, n.node, xi);
  apf::mapLocalToGlobal(me, xi, x);
  apf::destroyElement(e);
  apf::destroyMeshElement(me);
}

void StridedIndexer::compute_coords() {
  int dim = disc->get_num_dims();
  coords = rcp(new MultiVector(node_map, dim, false));
  apf::Vector3 point(0, 0, 0);
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(owned_numbering, owned);
  auto apf_f = fields[0]->get_apf_field();
  for (size_t n = 0; n < owned.size(); ++n) {
    get_coord(apf_f, owned[n], point);
    for (int d = 0; d < dim; ++d)
      coords->replaceLocalValue(n, d, point[d]);
  }
}

void StridedIndexer::compute_node_sets() {
  auto sets = disc->get_model_sets();
  auto nns = sets->models[0].size();
  auto nf = get_num_dof_fields();
  for (size_t s = 0; s < nns; ++s) {
    auto name = disc->get_node_set_name(s);
    node_sets[name].resize(nf);
    for (int f = 0; f < nf; ++f)
      node_sets[name][f].resize(0);
  }
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(owned_numbering, owned);
  for (size_t n = 0; n < owned.size(); ++n) {
    auto node = owned[n];
    auto ent = node.entity;
    std::set<apf::StkModel*> mset;
    apf::collectEntityModels(
        mesh, sets->invMaps[0], mesh->toModel(ent), mset);
    if (mset.empty()) continue;
    APF_ITERATE(std::set<apf::StkModel*>, mset, mit) {
        auto ns = *mit;
        auto nsn = ns->stkName;
        for (int f = 0; f < nf; ++f)
          node_sets[nsn][f].push_back(node);
    }
  }
}

} /* namespace goal */
