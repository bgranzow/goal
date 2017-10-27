#include "goal_assembly.hpp"
#include "goal_control.hpp"
#include "goal_dbcs.hpp"
#include "goal_disc.hpp"
#include "goal_eval_modes.hpp"
#include "goal_displacement.hpp"
#include "goal_linear_solve.hpp"
#include "goal_mechanics.hpp"
#include "goal_pressure.hpp"
#include "goal_primal.hpp"
#include "goal_scalar_weight.hpp"
#include "goal_sol_info.hpp"
#include "goal_vector_weight.hpp"

namespace goal {

using Teuchos::rcp;

static void make_displacement(Mechanics* m, Evaluators& r, Evaluators& j) {
  auto f = m->get_displacement();
  auto u = rcp(new Displacement<ST>(f, PRIMAL));
  auto up = rcp(new Displacement<FADT>(f, PRIMAL));
  auto w = rcp(new VectorWeight(f));
  r.push_back(u);
  r.push_back(w);
  j.push_back(up);
  j.push_back(w);
}

static void make_pressure(Mechanics* m, Evaluators& r, Evaluators& j) {
  auto f = m->get_pressure();
  auto p = rcp(new Pressure<ST>(f, PRIMAL));
  auto pp = rcp(new Pressure<FADT>(f, PRIMAL));
  auto w = rcp(new ScalarWeight(f));
  r.push_back(p);
  r.push_back(w);
  j.push_back(pp);
  j.push_back(w);
}

Primal::Primal(ParameterList const& p, Mechanics* m) {
  params = p;
  mech = m;
  sol_info = 0;
  make_displacement(mech, residual, jacobian);
  make_pressure(mech, residual, jacobian);
  mech->build_resid<ST>(residual, true);
  mech->build_resid<FADT>(jacobian, true);
}

Primal::~Primal() {
  destroy_data();
}

void Primal::build_data() {
  auto disc = mech->get_disc();
  sol_info = create_sol_info(disc);
}

void Primal::destroy_data() {
  if (sol_info) {
    destroy_sol_info(sol_info);
    sol_info = 0;
  }
}

void Primal::print_banner(double t_now) {
  auto ndofs = sol_info->owned->R->getGlobalLength();
  print("*** primal solve: %d dofs", ndofs);
  print("*** at time: %f", t_now);
}

void Primal::compute_resid(double t_now, double t_old) {
  auto t0 = time();
  auto dbc = params.sublist("dirichlet bcs");
  sol_info->zero_R();
  set_time(residual, t_now, t_old);
  assemble(residual, sol_info);
  sol_info->gather_R();
  set_resid_dbcs(dbc, sol_info, t_now);
  auto t1 = time();
  print(" > residual computed in %f seconds", t1 - t0);
}

void Primal::compute_jacob(double t_now, double t_old) {
  auto t0 = time();
  auto dbc = params.sublist("dirichlet bcs");
  sol_info->resume_fill();
  sol_info->zero_all();
  set_time(jacobian, t_now, t_old);
  assemble(jacobian, sol_info);
  sol_info->gather_all();
  set_jac_dbcs(dbc, sol_info, t_now);
  sol_info->complete_fill();
  auto t1 = time();
  print(" > jacobian computed in %f seconds", t1 - t0);
}

void Primal::solve(double t_now, double t_old) {
  print_banner(t_now);
  auto disc = sol_info->get_disc();
  auto R = sol_info->owned->R;
  auto dRdu = sol_info->owned->dRdu;
  auto du = rcp(new VectorT(disc->get_owned_map()));
  auto lp = params.sublist("primal linear algebra");
  auto max_iters = lp.get<int>("nonlinear max iters");
  auto tol = lp.get<double>("nonlinear tolerance");
  int iter = 1;
  bool converged = false;
  while ((iter <= max_iters) && (! converged)) {
    print(" > (%d) newton iteration", iter);
    compute_jacob(t_now, t_old);
    R->scale(-1.0);
    du->putScalar(0.0);
    goal::solve(lp, dRdu, du, R, disc);
    disc->add_soln(du);
    compute_resid(t_now, t_old);
    double norm = R->norm2();
    print(" > ||R|| = %e", norm);
    if (norm < tol) converged = true;
    iter++;
  }
  if ((iter > max_iters) && (! converged))
    fail("newton's method failed in %d iterations", max_iters);
}

Primal* create_primal(ParameterList const& p, Mechanics* m) {
  return new Primal(p, m);
}

void destroy_primal(Primal* p) {
  delete p;
}

}
