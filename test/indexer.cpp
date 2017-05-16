#include <apfMesh2.h>
#include <apfNumbering.h>
#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_field.hpp>
#include <goal_indexer.hpp>
#include <Teuchos_ParameterList.hpp>

namespace test {

static goal::Discretization* load_disc(char** argv) {
  Teuchos::ParameterList p;
  p.set<std::string>("geom file", argv[1]);
  p.set<std::string>("mesh file", argv[2]);
  p.set<std::string>("assoc file", argv[3]);
  p.set<bool>("reorder mesh", true);
  p.set<int>("workset size", 1000);
  return goal::create_disc(p);
}

static std::vector<goal::Field*> make_fields(goal::Discretization* d) {
  int dim = d->get_num_dims();
  std::vector<goal::Field*> u;
  if (dim > 0)
    u.push_back(goal::create_field({d, "ux", 1, 1, goal::LAGRANGE}));
  if (dim > 1)
    u.push_back(goal::create_field({d, "uy", 1, 1, goal::LAGRANGE}));
  if (dim > 2)
    u.push_back(goal::create_field({d, "uz", 1, 1, goal::LAGRANGE}));
  for (int i = 0; i < dim; ++i)
    u[i]->set_associated_dof_idx(i);
  return u;
}

static void destroy_fields(std::vector<goal::Field*>& u) {
  for (size_t i = 0; i < u.size(); ++i)
    goal::destroy_field(u[i]);
  u.resize(0);
}

static void check_tpetra_objs(goal::Indexer* i) {
  auto om = i->get_owned_map();
  auto gm = i->get_ghost_map();
  auto og = i->get_owned_graph();
  auto gg = i->get_ghost_graph();
  auto c = i->get_coords();
  goal::print("owned map elems: %lu", om->getGlobalNumElements());
  goal::print("ghost map elems: %lu", gm->getGlobalNumElements());
  goal::print("owned graph rows: %lu", og->getGlobalNumRows());
  goal::print("owned graph cols: %lu", og->getGlobalNumCols());
  goal::print("ghost graph rows: %lu", gg->getGlobalNumRows());
  goal::print("ghost graph cols: %lu", gg->getGlobalNumCols());
  goal::print("coords: num vecs: %lu", c->getNumVectors());
  goal::print("coords: num rows: %lu", c->getGlobalLength());
}

static void check_ids_by_type(goal::Indexer* i, apf::MeshEntity* e) {
  auto type = i->get_apf_mesh()->getType(e);
  static std::vector<goal::LO> lids;
  i->get_ghost_lids(e, lids);
  GOAL_ALWAYS_ASSERT(i->get_num_total_dofs(type) == (int)lids.size());
  for (size_t ldof = 0; ldof < lids.size(); ++ldof)
    GOAL_ALWAYS_ASSERT(lids[ldof] >= 0);
  for (int f = 0; f < i->get_num_fields(); ++f) {
    auto field = i->get_field(f);
    auto nnodes = field->get_num_nodes(type);
    for (int n = 0; n < nnodes; ++n)
      GOAL_ALWAYS_ASSERT(i->get_ghost_lid(f, e, n) >= 0);
  }
}

static void check_elem_ids(goal::Discretization* d, goal::Indexer* i) {
  goal::print("checking elem ids");
  for (int es = 0; es < d->get_num_elem_sets(); ++es) {
    auto name = d->get_elem_set_name(es);
    for (int ws = 0; ws < d->get_num_elem_worksets(es); ++ws) {
      auto elems = d->get_elems(name, ws);
      for (size_t e = 0; e < elems.size(); ++e)
        check_ids_by_type(i, elems[e]);
    }
  }
}

static void check_side_ids(goal::Discretization* d, goal::Indexer* i) {
  goal::print("checking side ids");
  for (int ss = 0; ss < d->get_num_side_sets(); ++ss) {
    auto name = d->get_side_set_name(ss);
    for (int ws = 0; ws < d->get_num_side_worksets(ss); ++ws) {
      auto sides = d->get_sides(name, ws);
      for (size_t s = 0; s < sides.size(); ++s)
        check_ids_by_type(i, sides[s]);
    }
  }
}

static void check_node_ids(goal::Discretization* d, goal::Indexer* i) {
  goal::print("checking node ids");
  for (int ns = 0; ns < d->get_num_node_sets(); ++ns) {
    auto name = d->get_node_set_name(ns);
    auto nodes = i->get_node_set_nodes(name);
    for (int f = 0; f < i->get_num_fields(); ++f) {
      for (size_t n = 0; n < nodes.size(); ++n)
        GOAL_ALWAYS_ASSERT(i->get_owned_lid(f, nodes[n]) >= 0);
    }
  }
}

void check_indexer(
    goal::Discretization* d, std::vector<goal::Field*> const& u) {
  goal::Indexer* i = goal::create_indexer(d, u);
  check_tpetra_objs(i);
  check_elem_ids(d, i);
  check_side_ids(d, i);
  check_node_ids(d, i);
  destroy_indexer(i);
}

} // end namespace test

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: strided indexer");
  GOAL_ALWAYS_ASSERT(argc == 4);
  auto d = test::load_disc(argv);
  auto u = test::make_fields(d);
  test::check_indexer(d, u);
  test::destroy_fields(u);
  goal::destroy_disc(d);
  goal::finalize();
}
