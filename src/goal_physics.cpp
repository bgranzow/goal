#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_assembly.hpp"
#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_physics.hpp"
#include "goal_point_wise.hpp"
#include "goal_residual.hpp"
#include "goal_scalar_types.hpp"

namespace goal {

using Teuchos::rcp;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("f", "");
  return p;
}

Physics::Physics(ParameterList const& p, Disc* d) {
  p.validateParameters(get_valid_params(), 0);
  params = p;
  disc = d;
  soln = 0;
  make_soln();
}

Physics::~Physics() {
  apf::destroyField(soln);
}

void Physics::make_soln() {
  auto m = disc->get_apf_mesh();
  auto f = m->findField("u");
  if (f) soln = f;
  else soln = apf::createFieldOn(m, "u", apf::SCALAR);
  if (!f) apf::zeroField(soln);
}

template <typename T>
void Physics::build_resid(Evaluators& E) {
  auto f = params.get<std::string>("f");
  auto u = find_evaluator("u", E);
  auto w = find_evaluator("uw", E);
  auto R = rcp(new Residual<T>(u, w, f));
  E.push_back(R);
}

template <typename T>
void Physics::build_functional(ParameterList const& params, Evaluators& E) {
  auto type = params.get<std::string>("type");
  auto u = find_evaluator("u", E);
  RCP<QoI<T>> J;
  if (type == "point wise")
    J = rcp(new PointWise<T>(params));
  else
    fail("unkown qoi: %s", type.c_str());
  E.push_back(J);
}

void Physics::build_error(Evaluators& E) {
  auto f = params.get<std::string>("f");
  auto u = find_evaluator("u", E);
  auto w = find_evaluator("uw", E);
  auto R = rcp(new Residual<ST>(u, w, f));
  E.push_back(R);
}

Physics* create_physics(ParameterList const& p, Disc* d) {
  return new Physics(p, d);
}

void destroy_physics(Physics* m) {
  delete m;
}

template void Physics::build_resid<ST>(Evaluators&);
template void Physics::build_resid<FADT>(Evaluators&);
template void Physics::build_functional<ST>(ParameterList const&, Evaluators&);
template void Physics::build_functional<FADT>(ParameterList const&, Evaluators&);

}
