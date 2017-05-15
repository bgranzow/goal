#include <apfAlbany.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_strided_indexer.hpp"

namespace goal {

static void validate_fields(std::vector<Field*> const& f) {
  for (size_t i = 0; i < f.size()-1; ++i) {
    auto b1 = f[i]->get_basis_type();
    auto b2 = f[i+1]->get_basis_type();
    auto p1 = f[i]->get_p_order();
    auto p2 = f[i+1]->get_p_order();
    GOAL_ALWAYS_ASSERT(b1 == b2);
    GOAL_ALWAYS_ASSERT(p1 == p2);
  }
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
    Discretization* d, std::vector<Field*> const& f)
    : Indexer(d, f),
      owned_nmbr(0),
      ghost_nmbr(0),
      global_nmbr(0) {
  auto t0 = time();
  validate_fields(fields);
  num_eqs = fields.size();
  mesh = disc->get_apf_mesh();
  shape = fields[0]->get_apf_basis();
  compute_owned_maps();
  compute_ghost_map();
  compute_graphs();
  compute_coords();
  compute_node_sets();
  auto t1 = time();
  print(" > strided indexer: pde equations: %d", num_eqs);
  print(" > strided indexer: built in %f seconds", t1 - t0);
}

StridedIndexer::~StridedIndexer() {
  destroy_local(&owned_nmbr);
  destroy_local(&ghost_nmbr);
  destroy_global(&global_nmbr);
}

static LO get_dof(const LO nid, const int eq, const int eqs) {
  return nid * eqs + eq;
}

static GO get_gdof(const LO nid, const int eq, const int eqs) {
  return nid * eqs + eq;
}

int StridedIndexer::get_elem_dof_offset(const int i, const int n) {
  return n * num_eqs + i;
}

LO StridedIndexer::get_ghost_lid(
    const int i, apf::MeshEntity* e, const int n) {
  static apf::NewArray<int> node_ids;
  apf::getElementNumbers(ghost_nmbr, e, node_ids);
  return get_dof(node_ids[n], i, num_eqs);
}

LO StridedIndexer::get_owned_lid(const int i, apf::Node const& n) {
  LO nid = apf::getNumber(owned_nmbr, n.entity, n.node, 0);
  return get_dof(nid, i, num_eqs);
}

void StridedIndexer::get_ghost_lids(
    apf::MeshEntity* e, std::vector<LO>& lids) {
  auto t = mesh->getType(e);
  lids.resize(get_num_total_dofs(t));
  static apf::NewArray<int> nids;
  apf::getElementNumbers(ghost_nmbr, e, nids);
  int dof = 0;
  for (int n = 0; n < fields[0]->get_num_nodes(t); ++n)
    for (int eq = 0; eq < num_eqs; ++eq)
      lids[dof++] = get_dof(nids[n], eq, num_eqs);
}

void StridedIndexer::add_to_fields(
    std::vector<Field*> const& flds, RCP<Vector> du) {
  (void)flds;
  (void)du;
}

void StridedIndexer::set_to_fields(
    std::vector<Field*> const& flds, RCP<Vector> x) {
  (void)flds;
  (void)x;
}

void StridedIndexer::compute_owned_maps() {
  destroy_local(&owned_nmbr);
  destroy_global(&global_nmbr);
  owned_nmbr = apf::numberOwnedNodes(mesh, "owned", shape);
  global_nmbr = apf::makeGlobal(owned_nmbr, false);
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(global_nmbr, owned);
  auto num_owned = owned.getSize();
  Teuchos::Array<GO> indices(num_owned);
  for (size_t n = 0; n < num_owned; ++n)
    indices[n] = apf::getNumber(global_nmbr, owned[n]);
  node_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
  indices.resize(num_eqs * num_owned);
  for (size_t n = 0; n < num_owned; ++n) {
    GO gid = apf::getNumber(global_nmbr, owned[n]);
    for (int eq = 0; eq < num_eqs; ++eq)
      indices[get_dof(n, eq, num_eqs)] = get_gdof(gid, eq, num_eqs);
  }
  owned_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
  apf::synchronize(global_nmbr);
}

void StridedIndexer::compute_ghost_map() {
  destroy_local(&ghost_nmbr);
  ghost_nmbr = apf::numberOverlapNodes(mesh, "ghost", shape);
  apf::DynamicArray<apf::Node> ghost;
  apf::getNodes(ghost_nmbr, ghost);
  auto num_ghost = ghost.getSize();
  Teuchos::Array<GO> indices(num_eqs * num_ghost);
  for (size_t n = 0; n < num_ghost; ++n) {
    GO gid = apf::getNumber(global_nmbr, ghost[n]);
    for (int eq = 0; eq < num_eqs; ++eq)
      indices[get_dof(n, eq, num_eqs)] = get_gdof(gid, eq, num_eqs);
  }
  ghost_map = Tpetra::createNonContigMap<LO, GO>(indices, comm);
}

void StridedIndexer::compute_graphs() {
  int dim = disc->get_num_dims();
  int est = get_num_total_dofs(apf::Mesh::TET);
  owned_graph = rcp(new Graph(owned_map, est));
  ghost_graph = rcp(new Graph(ghost_map, est));
  apf::MeshEntity* elem;
  auto elems = mesh->begin(dim);
  while ((elem = mesh->iterate(elems))) {
    apf::NewArray<long> gids;
    int nnodes = apf::getElementNumbers(global_nmbr, elem, gids);
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
  destroy_global(&global_nmbr);
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
  apf::getNodes(owned_nmbr, owned);
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
  auto nf = get_num_fields();
  for (size_t s = 0; s < nns; ++s) {
    auto name = disc->get_node_set_name(s);
    node_sets[name].resize(nf);
    for (int f = 0; f < nf; ++f)
      node_sets[name][f].resize(nf);
  }
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(owned_nmbr, owned);
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

} // end namespace goal
