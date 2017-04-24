#include <apfAlbany.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

#include "goal_control.hpp"
#include "goal_field.hpp"
#include "goal_strided_indexer.hpp"

namespace goal {

static int get_num_eqs(std::vector<RCP<Field> > const& f) {
  int num_eqs = 0;
  for (size_t i = 0; i < f.size(); ++i)
    num_eqs += f[i]->get_num_components();
  return num_eqs;
}

StridedIndexer::StridedIndexer(
    RCP<Discretization> d, std::vector<RCP<Field> > f)
    : Indexer(d, f) {
  auto t0 = time();
  num_eqs = get_num_eqs(fields);
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
  if (owned_numbering) apf::destroyNumbering(owned_numbering);
  if (ghost_numbering) apf::destroyNumbering(owned_numbering);
  if (global_numbering) apf::destroyGlobalNumbering(global_numbering);
}

static LO get_dof(const LO node_id, const LO eq, const LO num_eqs) {
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
    std::vector<RCP<Field> > const& f, RCP<Vector> du) {
  (void)f;
  (void)du;
}

void StridedIndexer::set_to_fields(
    std::vector<RCP<Field> > const& f, RCP<Vector> x) {
  (void)f;
  (void)x;
}

void StridedIndexer::compute_owned_maps() {
}

void StridedIndexer::compute_ghost_map() {
}

void StridedIndexer::compute_graphs() {
}

void StridedIndexer::compute_coords() {
}

void StridedIndexer::compute_node_sets() {
}

} /* namespace goal */
