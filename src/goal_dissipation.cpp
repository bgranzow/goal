#include "goal_assembly.hpp"
#include "goal_control.hpp"
#include "goal_dbcs.hpp"
#include "goal_disc.hpp"
#include "goal_displacement.hpp"
#include "goal_dissipation.hpp"
#include "goal_eval_modes.hpp"
#include "goal_linear_solve.hpp"
#include "goal_mechanics.hpp"
#include "goal_pressure.hpp"
#include "goal_primal.hpp"
#include "goal_scalar_weight.hpp"
#include "goal_sol_info.hpp"
#include "goal_vector_weight.hpp"

#include <apf.h>
#include <apfMesh2.h>

namespace goal {

using Teuchos::rcp;

static void make_displacement_adj(Mechanics* m, Evaluators& a) {
  auto f = m->get_displacement();
  auto u = rcp(new Displacement<FADT>(f, ADJOINT));
  auto w = rcp(new VectorWeight(f));
  a.push_back(u);
  a.push_back(w);
}

static void make_pressure_adj(Mechanics* m, Evaluators& a) {
  auto f = m->get_pressure();
  auto p = rcp(new Pressure<FADT>(f, ADJOINT));
  auto w = rcp(new ScalarWeight(f));
  a.push_back(p);
  a.push_back(w);
}

Dissipation::Dissipation(ParameterList const& p, Primal* pr) {
  params = p;
  primal = pr;
  mech = 0;
  sol_info = 0;
  disc = 0;
}

Dissipation::~Dissipation() {
}

void Dissipation::print_banner(double t_now) {
  print("*** dissipation error indication");
  print("*** at time: %f", t_now);
}

void Dissipation::build_adjoint() {
  auto func_params = params.sublist("functional");
  sol_info = primal->get_sol_info();
  mech = primal->get_mech();
  disc = mech->get_disc();
  make_displacement_adj(mech, adjoint);
  make_pressure_adj(mech, adjoint);
  mech->build_resid<FADT>(adjoint, false);
  mech->build_functional<FADT>(func_params, adjoint);
  auto mesh = disc->get_apf_mesh();
  z_displacement = apf::createFieldOn(mesh, "zu", apf::VECTOR);
  z_pressure = apf::createFieldOn(mesh, "zp", apf::SCALAR);
}

void Dissipation::destroy_adjoint() {
  apf::destroyField(z_displacement);
  apf::destroyField(z_pressure);
  z_displacement = 0;
  z_pressure = 0;
  sol_info = 0;
  mech = 0;
  disc = 0;
  adjoint.resize(0);
}

void Dissipation::compute_adjoint(double t_now, double t_old) {
  auto t0 = time();
  auto dbc = params.sublist("dirichlet bcs");
  sol_info->resume_fill();
  sol_info->zero_all();
  set_time(adjoint, t_now, t_old);
  assemble(adjoint, sol_info);
  sol_info->gather_all();
  set_jac_dbcs(dbc, sol_info, t_now);
  sol_info->complete_fill();
  auto t1 = time();
  print(" > adjoint computed in %f seconds", t1 - t0);
}

void Dissipation::solve_adjoint(double t_now, double t_old) {
  print_banner(t_now);
  auto dRduT = sol_info->owned->dRdu;
  auto dMdu = sol_info->owned->dMdu;
  auto z = rcp(new VectorT(disc->get_owned_map()));
  auto lp = params.sublist("adjoint linear algebra");
  compute_adjoint(t_now, t_old);
  z->putScalar(0.0);
  goal::solve(lp, dRduT, z, dMdu, disc);
  disc->set_adj(z);
}

void Dissipation::run(double t_now, double t_old) {
  build_adjoint();
  solve_adjoint(t_now, t_old);
  apf::writeVtkFiles("DEBUG",  disc->get_apf_mesh());
  destroy_adjoint();
}

Dissipation* create_dissipation(ParameterList const& p, Primal* pr) {
  return new Dissipation(p, pr);
}

void destroy_dissipation(Dissipation* d) {
  delete d;
}

}
