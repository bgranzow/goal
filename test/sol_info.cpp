#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_field.hpp>
#include <goal_indexer.hpp>
#include <goal_sol_info.hpp>
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

static void check_sol_info(
    goal::Discretization* d, std::vector<goal::Field*> const& u) {
  goal::Indexer* i = goal::create_indexer(d, u);
  goal::SolInfo* s = goal::create_sol_info(i);
  goal::print("owned du rows: %lu", s->owned->du->getGlobalLength());
  goal::print("ghost du rows: %lu", s->ghost->du->getGlobalLength());
  goal::print("owned R rows: %lu", s->owned->R->getGlobalLength());
  goal::print("ghost R rows: %lu", s->ghost->R->getGlobalLength());
  goal::print("owned dRdu rows: %lu", s->owned->dRdu->getGlobalNumRows());
  goal::print("ghost dRdu cols: %lu", s->ghost->dRdu->getGlobalNumCols());
  goal::print("owned dJdu rows: %lu", s->owned->dJdu->getGlobalLength());
  goal::print("ghost dJdu rows: %lu", s->ghost->dJdu->getGlobalLength());
  goal::destroy_sol_info(s);
  goal::destroy_indexer(i);
}

}

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: solution info");
  GOAL_ALWAYS_ASSERT(argc == 4);
  auto d = test::load_disc(argv);
  auto u = test::make_fields(d);
  test::check_sol_info(d, u);
  test::destroy_fields(u);
  goal::destroy_disc(d);
  goal::finalize();
}
