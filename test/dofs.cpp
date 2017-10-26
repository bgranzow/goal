#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include <goal_control.hpp>
#include <goal_disc.hpp>
#include <goal_displacement.hpp>
#include <goal_eval_modes.hpp>
#include <goal_pressure.hpp>
#include <goal_sol_info.hpp>
#include <Teuchos_ParameterList.hpp>

namespace test {

using Teuchos::rcp;

template <typename T>
void eval(goal::SolInfo* s, apf::Field* uf, apf::Field* pf) {
  apf::Vector3 point;
  auto u = rcp(new goal::Displacement<T>(uf, goal::PRIMAL));
  auto p = rcp(new goal::Pressure<T>(pf, goal::PRIMAL));
  auto d = s->get_disc();
  auto m = d->get_apf_mesh();
  u->pre_process(s);
  p->pre_process(s);
  for (int es = 0; es < d->get_num_elem_sets(); ++es) {
    auto esn = d->get_elem_set_name(es);
    auto elems = d->get_elems(esn);
    for (size_t elem = 0; elem < elems.size(); ++elem) {
      auto me = apf::createMeshElement(m, elems[elem]);
      u->gather(me);
      p->gather(me);
      apf::getIntPoint(me, 1, 0, point);
      u->at_point(point, 0, 0);
      p->at_point(point, 0, 0);
      u->out_elem();
      p->out_elem();
      apf::destroyMeshElement(me);
    }
  }
  u->post_process(NULL);
  p->post_process(NULL);
}

static void check_dofs(goal::SolInfo* s) {
  auto m = s->get_disc()->get_apf_mesh();
  auto u = apf::createFieldOn(m, "u", apf::VECTOR);
  auto p = apf::createFieldOn(m, "p", apf::VECTOR);
  apf::zeroField(u);
  apf::zeroField(p);
  eval<goal::ST>(s, u, p);
  eval<goal::FADT>(s, u, p);
  apf::destroyField(p);
  apf::destroyField(u);
}

}

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: dofs");
  GOAL_ALWAYS_ASSERT(argc == 4);
  Teuchos::ParameterList p;
  p.set<std::string>("geom file", argv[1]);
  p.set<std::string>("mesh file", argv[2]);
  p.set<std::string>("assoc file", argv[3]);
  auto d = goal::create_disc(p);
  d->build_data();
  auto s = goal::create_sol_info(d);
  test::check_dofs(s);
  goal::destroy_sol_info(s);
  goal::destroy_disc(d);
  goal::finalize();
}
