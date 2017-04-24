#include <apfAlbany.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

#include "goal_control.hpp"
#include "goal_field.hpp"
#include "goal_strided_indexer.hpp"

namespace goal {

static void get_num_eqs(std::vector<RCP<Field> > const& f) {
  int num_eqs = 0;
  for (size_t i = 0; i < f.size(); ++i)
    num_eqs += f->get_num_components();
  return num_eqs;
}

StridedIndexer::StridedIndexer(
    RCP<Discretization> d, std::vector<RCP<Field> > f)
    : Indexer(d, f) {
  auto t0 = time();
  num_eqs = get_num_eqs(fields);
  compute_owned_map();
  compute_ghost_map();
  compute_graphs();
  compute_node_map();
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

int StridedIndexer::get_elem_dof_offset(
    const int idx, const int n, const int c) {
  int eq = 0;
  for (int i = 0; i < idx; ++i)
    eq += fields[i]->get_num_components();
  eq += c;
  return n*num_eqs + eq;
}

LO StridedIndexer::get_ghost_lid(
    const int idx, apf::MeshEntity* e, const int n, const int c) {
}

LO StridedIndexer::get_owned_lid(
    const int idx, apf::Node const& n, const int c) {
}

void StridedIndexer::get_ghost_lids(
    apf::MeshEntity* e, std::vector<LO>& lids) {
}

void StridedIndexer::add_to_fields(
    std::vector<RCP<Field> > const& f, RCP<Vector> du) {
}

void StridedIndexer::set_to_fields(
    std::vector<RCP<Field> > const& f, RCP<Vector> x) {
}

void StridedIndexer::compute_owned_maps() {
}

void StridedIndexer::compute_ghost_map() {
}

void StridedIndexer::compute_graphs() {
}

void StridedIndexer:;compute_coords() {
}

void StridedIndexer::compute_node_sets() {
}

};
