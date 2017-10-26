#include <apfNumbering.h>
#include <goal_control.hpp>
#include <goal_disc.hpp>
#include <goal_nested.hpp>
#include <Teuchos_ParameterList.hpp>

namespace test {

static void check_sets(goal::Disc* d) {
  goal::print("num dims: %d", d->get_num_dims());
  goal::print("num elem sets: %d", d->get_num_elem_sets());
  goal::print("num side sets: %d", d->get_num_side_sets());
  goal::print("num node sets: %d", d->get_num_node_sets());
  for (int i = 0; i < d->get_num_elem_sets(); ++i) {
    goal::print("elem set: %d", i);
    goal::print(" name: %s", (d->get_elem_set_name(i)).c_str());
    d->get_elems(d->get_elem_set_name(i));
  }
  for (int i = 0; i < d->get_num_side_sets(); ++i) {
    goal::print("side set: %d", i);
    goal::print(" name: %s", (d->get_side_set_name(i)).c_str());
    d->get_sides(d->get_side_set_name(i));
  }
  for (int i = 0; i < d->get_num_node_sets(); ++i) {
    goal::print("node set: %d", i);
    goal::print(" name: %s", (d->get_node_set_name(i)).c_str());
    d->get_nodes(d->get_node_set_name(i));
  }
}

static void check_tpetra_objs(goal::Disc* d) {
  auto om = d->get_owned_map();
  auto gm = d->get_ghost_map();
  auto og = d->get_owned_graph();
  auto gg = d->get_ghost_graph();
  auto c = d->get_coords();
  goal::print("owned map size: %lu", om->getGlobalNumElements());
  goal::print("ghost map size: %lu", gm->getGlobalNumElements());
  goal::print("coords num vecs: %lu", c->getNumVectors());
  goal::print("coords num rows: %lu", c->getGlobalLength());
  goal::print("owned graph rows: %lu", og->getGlobalNumRows());
  goal::print("owned graph cols: %lu", og->getGlobalNumCols());
  goal::print("ghost graph rows: %lu", gg->getGlobalNumRows());
  goal::print("ghost graph cols: %lu", gg->getGlobalNumCols());
}

static void check_elem_indices(goal::Disc* d) {
  auto es_name = d->get_elem_set_name(0);
  auto elem = d->get_elems(es_name)[0];
  auto num_nodes = d->get_num_nodes(elem);
  auto num_dofs = d->get_num_dofs(elem);
  auto num_eqs = d->get_num_eqs();
  goal::print("num nodes: %d", num_nodes);
  goal::print("num dofs: %d", num_dofs);
  for (int n = 0; n < num_nodes; ++n)
  for (int eq = 0; eq < num_eqs; ++eq)
    goal::print("lid: %d", d->get_lid(elem, n, eq));
  std::vector<goal::LO> lids;
  d->get_lids(elem, lids);
  for (int i = 0; i < num_dofs; ++i)
    goal::print("lids[%d]: %d", i, lids[i]);
}

static void check_node_indices(goal::Disc* d) {
  auto ns_name = d->get_node_set_name(0);
  auto nodes = d->get_nodes(ns_name);
  if (nodes.size() == 0) return;
  int num_eqs = d->get_num_eqs();
  for (int eq = 0; eq < num_eqs; ++eq)
    goal::print("node id: %d", d->get_lid(nodes[0], eq));
}

static void check_rebuild(goal::Disc* d) {
  d->destroy_data();
  d->build_data();
}

static void check_disc(goal::Disc* d) {
  d->build_data();
  check_sets(d);
  check_tpetra_objs(d);
  check_elem_indices(d);
  check_node_indices(d);
  check_rebuild(d);
}

}

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: disc");
  GOAL_ALWAYS_ASSERT(argc == 4);
  Teuchos::ParameterList p;
  p.set<std::string>("geom file", argv[1]);
  p.set<std::string>("mesh file", argv[2]);
  p.set<std::string>("assoc file", argv[3]);
  auto d = goal::create_disc(p);
  auto n1 = goal::create_nested(d, 0);
  if (0)
  test::check_disc(d);
  goal::destroy_nested(n1);
  goal::destroy_disc(d);
  goal::finalize();
}
