#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_field.hpp>
#include <goal_output.hpp>
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

static goal::Field* make_field(goal::Discretization* d) {
  goal::FieldInfo i = {d, "u", 2, 2, goal::LAGRANGE}; 
  auto f = goal::create_field(i);
  f->set_seed(1.0);
  f->set_associated_dof_idx(0);
  return f;
}

static goal::Output* make_output(goal::Discretization* d) {
  Teuchos::ParameterList p;
  Teuchos::Array<std::string> f(1);
  f[0] = "u";
  p.set<std::string>("out file", "test_out");
  p.set<int>("interval", 2);
  p.set<Teuchos::Array<std::string> >("interpolate", f);
  auto out = goal::create_output(p, d);
  return out;
}

} // end namespace test

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: output");
  GOAL_ALWAYS_ASSERT(argc == 4);
  Teuchos::ParameterList p;
  p.set<std::string>("geom file", argv[1]);
  p.set<std::string>("mesh file", argv[2]);
  p.set<std::string>("assoc file", argv[3]);
  p.set<bool>("reorder mesh", true);
  p.set<int>("workset size", 1000);
  auto d = test::load_disc(argv);
  auto f = test::make_field(d);
  auto o = test::make_output(d);
  for (int i = 0; i < 4; ++i)
    o->write((double)i);
  goal::destroy_output(o);
  goal::destroy_field(f);
  goal::destroy_disc(d);
  goal::finalize();
}
