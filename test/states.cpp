#include <apf.h>
#include <apfMesh2.h>
#include <goal_control.hpp>
#include <goal_data_types.hpp>
#include <goal_discretization.hpp>
#include <goal_states.hpp>
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

static void check_bulk(goal::States* s, goal::Output* o) {
  s->add("cauchy", 2);
  s->add("eqps", 0, true);
  s->add("Fp", 2, true, true);
  o->write(0);
  s->project(2);
  o->write(1);
  s->project(1);
  o->write(2);
  s->update();
  o->write(3);
}

template <typename T>
void check_ent(goal::Discretization* d, goal::States* s) {
  using Tensor = minitensor::Tensor<T>;
  auto m = d->get_apf_mesh();
  auto dim = d->get_num_dims();
  apf::MeshEntity* e;
  apf::MeshIterator* elems = m->begin(dim);
  T eqps;
  Tensor Fp(dim);
  Tensor sigma(dim);
  Tensor I(minitensor::eye<T>(dim));
  while ((e = m->iterate(elems))) {
    s->get_scalar("eqps", e, 0, eqps);
    s->get_tensor("Fp", e, 0, Fp);
    s->get_tensor("cauchy", e, 0, sigma);
    eqps += 1.0;
    Fp += I;
    sigma += I;
    s->set_scalar("eqps", e, 0, eqps);
    s->set_tensor("Fp", e, 0, Fp);
    s->set_tensor("cauchy", e, 0, sigma);
  }
  m->end(elems);
}

static void check_states(goal::Discretization* d) {
  Teuchos::ParameterList p;
  p.set<std::string>("out file", "states");
  auto o = goal::create_output(p, d);
  auto s = goal::create_states(d, 1);
  check_bulk(s, o);
  check_ent<double>(d, s);
  check_ent<goal::FadType>(d, s);
  o->write(4);
  s->update();
  o->write(5);
  goal::destroy_states(s);
  goal::destroy_output(o);
}

} // end namespace test

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: states");
  GOAL_ALWAYS_ASSERT(argc == 4);
  auto d = test::load_disc(argv);
  test::check_states(d);
  goal::destroy_disc(d);
  goal::finalize();
}
