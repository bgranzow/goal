#include <apf.h>
#include <apfMesh2.h>

#include <goal_control.hpp>
#include <goal_disc.hpp>
#include <goal_scalar_types.hpp>
#include <goal_states.hpp>

namespace test {

static void check_bulk(goal::States* s) {
  s->add("cauchy", apf::MATRIX);
  s->add("eqps", apf::SCALAR, true);
  s->add("Fp", apf::MATRIX, true, true);
  s->update();
}

template <typename T>
static void check_ent(goal::States* s) {
  using Tensor = minitensor::Tensor<T>;
  auto d = s->get_disc();
  auto m = d->get_apf_mesh();
  auto dim = d->get_num_dims();
  T eqps;
  Tensor Fp(dim);
  Tensor sigma(dim);
  Tensor I(minitensor::eye<T>(dim));
  apf::MeshEntity* e;
  apf::MeshIterator* elems = m->begin(dim);
  while ((e = m->iterate(elems))) {
    goal::get_scalar(d, "eqps_old", e, eqps);
    goal::get_tensor(d, "Fp_old", e, Fp);
    goal::get_tensor(d, "cauchy", e, sigma);
    eqps += 1.0;
    Fp += I;
    sigma += I;
    goal::set_scalar(d, "eqps", e, eqps);
    goal::set_tensor(d, "Fp", e, Fp);
    goal::set_tensor(d, "cauchy", e, sigma);
  }
  m->end(elems);
}

static void check_states(goal::States* s) {
  check_bulk(s);
  check_ent<goal::ST>(s);
  check_ent<goal::FADT>(s);
}

}

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: states");
  GOAL_ALWAYS_ASSERT(argc == 4);
  Teuchos::ParameterList p;
  p.set<std::string>("geom file", argv[1]);
  p.set<std::string>("mesh file", argv[2]);
  p.set<std::string>("assoc file", argv[3]);
  auto d = goal::create_disc(p);
  auto s = goal::create_states(d);
  test::check_states(s);
  goal::destroy_states(s);
  goal::destroy_disc(d);
  goal::finalize();
}
