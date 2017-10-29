#include "goal_assembly.hpp"
#include "goal_control.hpp"
#include "goal_dbcs.hpp"
#include "goal_displacement.hpp"
#include "goal_eval_modes.hpp"
#include "goal_linear_solve.hpp"
#include "goal_mechanics.hpp"
#include "goal_nested.hpp"
#include "goal_nested_adjoint.hpp"
#include "goal_pressure.hpp"
#include "goal_primal.hpp"
#include "goal_scalar_weight.hpp"
#include "goal_sol_info.hpp"
#include "goal_tbcs.hpp"
#include "goal_vector_weight.hpp"

#include <apf.h>
#include <apfMesh2.h>

namespace goal {

using Teuchos::rcp;

static void make_displacement_adj(Disc* d, Evaluators& a) {
  auto m = d->get_apf_mesh();
  auto f = m->findField("u");
  GOAL_DEBUG_ASSERT(f);
  auto u = rcp(new Displacement<FADT>(f, ADJOINT));
  auto w = rcp(new VectorWeight(f));
  a.push_back(u);
  a.push_back(w);
}

static void make_pressure_adj(Disc* d, Evaluators& a) {
  auto m = d->get_apf_mesh();
  auto f = m->findField("p");
  GOAL_DEBUG_ASSERT(f);
  auto p = rcp(new Pressure<FADT>(f, ADJOINT));
  auto w = rcp(new ScalarWeight(f));
  a.push_back(p);
  a.push_back(w);
}

NestedAdjoint::NestedAdjoint(ParameterList const& p, Primal* pr) {
  params = p;
  primal = pr;
  base_disc = primal->get_mech()->get_disc();
  nested_disc = 0;
  mech = 0;
  sol_info = 0;
}

NestedAdjoint::~NestedAdjoint() {
  destroy_data();
}

static int get_mode(std::string const& m) {
  if (m == "uniform") return UNIFORM;
  else if (m == "long") return LONG;
  else fail("unknown adjoint mode: %s", m.c_str());
}

void NestedAdjoint::build_data() {
  auto mech_params = params.sublist("mechanics");
  auto func_params = params.sublist("functional");
  auto nested_mode = params.get<std::string>("adjoint mode");
  auto mode = get_mode(nested_mode);
  nested_disc = create_nested(base_disc, mode);
  mech = create_mechanics(mech_params, nested_disc);
  nested_disc->build_data();
  sol_info = create_sol_info(nested_disc);
  make_displacement_adj(nested_disc, adjoint);
  make_pressure_adj(nested_disc, adjoint);
  mech->build_resid<FADT>(adjoint, false);
  mech->build_functional<FADT>(func_params, adjoint);
}

void NestedAdjoint::destroy_data() {
  if (sol_info) destroy_sol_info(sol_info);
  if (mech) destroy_mechanics(mech);
  if (nested_disc) destroy_nested(nested_disc);
  nested_disc = 0;
  mech = 0;
  sol_info = 0;
}

void NestedAdjoint::print_banner(const double t_now) {
  auto ndofs = sol_info->owned->R->getGlobalLength();
  print("*** adjoint error estimation: %d dofs", ndofs);
  print("*** at time: %f", t_now);
}

void NestedAdjoint::compute_adjoint(double t_now, double t_old) {
  auto t0 = time();
  auto dbc = params.sublist("dirichlet bcs");
  auto tbc = params.sublist("traction bcs");
  auto w = find_evaluator("uw", adjoint);
  sol_info->resume_fill();
  sol_info->zero_all();
  set_time(adjoint, t_now, t_old);
  assemble(adjoint, sol_info);
  set_tbcs(tbc, w, sol_info, t_now);
  sol_info->gather_all();
  set_jac_dbcs(dbc, sol_info, t_now);
  sol_info->complete_fill();
  auto t1 = time();
  print(" > adjoint computed in %f seconds", t1 - t0);
}

void NestedAdjoint::solve(double t_now, double t_old) {
  auto R = sol_info->owned->R;
  auto dRduT = sol_info->owned->dRdu;
  auto dMdu = sol_info->owned->dMdu;
  auto z = rcp(new VectorT(nested_disc->get_owned_map()));
  auto lp = params.sublist("adjoint linear algebra");
  compute_adjoint(t_now, t_old);
  z->putScalar(0.0);
  goal::solve(lp, dRduT, z, dMdu, nested_disc);
}

void NestedAdjoint::run(double t_now, double t_old) {
  print_banner(t_now);
  solve(t_now, t_old);
}

NestedAdjoint* create_nested_adjoint(ParameterList const& p, Primal* pr) {
  return new NestedAdjoint(p, pr);
}

void destroy_nested_adjoint(NestedAdjoint* a) {
  delete a;
}

}
