#include "goal_assembly.hpp"
#include "goal_control.hpp"
#include "goal_dbcs.hpp"
#include "goal_displacement.hpp"
#include "goal_displacement_adjoint.hpp"
#include "goal_error.hpp"
#include "goal_eval_modes.hpp"
#include "goal_ibcs.hpp"
#include "goal_linear_solve.hpp"
#include "goal_mechanics.hpp"
#include "goal_nested.hpp"
#include "goal_nested_adjoint.hpp"
#include "goal_pressure.hpp"
#include "goal_pressure_adjoint.hpp"
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

static void make_displacement_err(Disc* d, Evaluators& e) {
  auto m = d->get_apf_mesh();
  auto fu = m->findField("u");
  auto fz = m->findField("u_z_diff");
  GOAL_DEBUG_ASSERT(fu);
  GOAL_DEBUG_ASSERT(fz);
  auto u = rcp(new Displacement<ST>(fu, PRIMAL));
  auto w = rcp(new DisplacementAdjoint(fz));
  e.push_back(u);
  e.push_back(w);
}

static void make_pressure_err(Disc* d, Evaluators& e) {
  auto m = d->get_apf_mesh();
  auto fp = m->findField("p");
  auto fz = m->findField("p_z_diff");
  auto fzc = m->findField("p_z_coarse");
  GOAL_DEBUG_ASSERT(fp);
  GOAL_DEBUG_ASSERT(fz);
  auto p = rcp(new Pressure<ST>(fp, PRIMAL));
  auto w = rcp(new PressureAdjoint(fz, "pw"));
  auto wc = rcp(new PressureAdjoint(fzc, "pwc"));
  e.push_back(p);
  e.push_back(w);
  e.push_back(wc);
}

NestedAdjoint::NestedAdjoint(ParameterList const& p, Primal* pr) {
  params = p;
  primal = pr;
  base_disc = primal->get_mech()->get_disc();
  nested_disc = 0;
  mech = 0;
  sol_info = 0;
  zu_coarse = 0;
  zu_fine = 0;
  zu_diff = 0;
  u_error = 0;
  zp_coarse = 0;
  zp_fine = 0;
  zp_diff = 0;
  p_error = 0;
}

NestedAdjoint::~NestedAdjoint() {
  destroy_data();
}

static int get_mode(std::string const& m) {
  if (m == "uniform") return UNIFORM;
  else if (m == "long") return LONG;
  else if (m == "single") return SINGLE;
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
  auto nested_mesh = nested_disc->get_apf_mesh();
  zu_coarse = apf::createFieldOn(nested_mesh, "u_z_coarse", apf::VECTOR);
  zu_fine = apf::createFieldOn(nested_mesh, "u_z_fine", apf::VECTOR);
  zu_diff = apf::createFieldOn(nested_mesh, "u_z_diff", apf::VECTOR);
  u_error = apf::createFieldOn(nested_mesh, "u_error", apf::VECTOR);
  zp_coarse = apf::createFieldOn(nested_mesh, "p_z_coarse", apf::SCALAR);
  zp_fine = apf::createFieldOn(nested_mesh, "p_z_fine", apf::SCALAR);
  zp_diff = apf::createFieldOn(nested_mesh, "p_z_diff", apf::SCALAR);
  p_error = apf::createFieldOn(nested_mesh, "p_error", apf::SCALAR);
  make_displacement_adj(nested_disc, adjoint);
  make_pressure_adj(nested_disc, adjoint);
  mech->build_resid<FADT>(adjoint, false);
  mech->build_functional<FADT>(func_params, adjoint);
  make_displacement_err(nested_disc, error);
  make_pressure_err(nested_disc, error);
  mech->build_error(error);
}

void NestedAdjoint::destroy_data() {
  if (p_error) apf::destroyField(p_error);
  if (zp_diff) apf::destroyField(zp_diff);
  if (zp_fine) apf::destroyField(zp_fine);
  if (zp_coarse) apf::destroyField(zp_coarse);
  if (u_error) apf::destroyField(u_error);
  if (zu_diff) apf::destroyField(zu_diff);
  if (zu_fine) apf::destroyField(zu_fine);
  if (zu_coarse) apf::destroyField(zu_coarse);
  if (sol_info) destroy_sol_info(sol_info);
  if (mech) destroy_mechanics(mech);
  if (nested_disc) destroy_nested(nested_disc);
  nested_disc = 0;
  mech = 0;
  sol_info = 0;
  zu_coarse = 0;
  zu_fine = 0;
  zu_diff = 0;
  u_error = 0;
  zp_coarse = 0;
  zp_fine = 0;
  zp_diff = 0;
  p_error = 0;
  adjoint.resize(0);
  error.resize(0);
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
  auto ibc = params.sublist("inward bcs");
  auto w = find_evaluator("uw", adjoint);
  sol_info->resume_fill();
  sol_info->zero_all();
  set_time(adjoint, t_now, t_old);
  assemble(adjoint, sol_info);
  set_tbcs(tbc, w, sol_info, t_now);
  set_ibcs(ibc, w, sol_info);
  sol_info->gather_all();
  set_jac_dbcs(dbc, sol_info, t_now);
  sol_info->complete_fill();
  auto t1 = time();
  print(" > adjoint computed in %f seconds", t1 - t0);
}

void NestedAdjoint::subtract() {
  apf::Vector3 zuf(0,0,0);
  apf::Vector3 zuc(0,0,0);
  auto m = nested_disc->get_apf_mesh();
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = m->begin(0);
  while ((vtx = m->iterate(it))) {
    apf::getVector(zu_fine, vtx, 0, zuf);
    apf::getVector(zu_coarse, vtx, 0, zuc);
    auto zpf = apf::getScalar(zp_fine, vtx, 0);
    auto zpc = apf::getScalar(zp_coarse, vtx, 0);
    apf::setVector(zu_diff, vtx, 0, zuf - zuc);
    apf::setScalar(zp_diff, vtx, 0, zpf - zpc);
  }
  m->end(it);
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
  nested_disc->set_fine(z, zu_fine, zp_fine);
  apf::copyData(zu_coarse, zu_fine);
  apf::copyData(zp_coarse, zp_fine);
  nested_disc->set_coarse(zu_coarse, zp_coarse);
  subtract();
  auto e = -(R->dot(*z));
  print(" > z.R = %.15e", e);
}

void NestedAdjoint::localize(double t_now, double t_old) {
  auto t0 = time();
  auto dbc = params.sublist("dirichlet bcs");
  auto tbc = params.sublist("traction bcs");
  auto ibc = params.sublist("inward bcs");
  auto w = find_evaluator("uw", error);
  auto R = sol_info->owned->R;
  sol_info->zero_R();
  set_time(error, t_now, t_old);
  assemble(error, sol_info);
  set_tbcs(tbc, w, sol_info, t_now);
  set_ibcs(ibc, w, sol_info);
  sol_info->gather_all();
  set_resid_dbcs(dbc, sol_info, t_now);
  nested_disc->set_fine(R, u_error, p_error);
  auto t1 = time();
  print(" > error localized in %f seconds", t1 - t0);
}

void NestedAdjoint::write_out() {
  static int i = 0;
  auto out_params = params.sublist("output");
  auto name = out_params.get<std::string>("out file");
  auto m = nested_disc->get_apf_mesh();
  std::ostringstream oss;
  oss << name << "_adjoint_" << i;
  std::string fname = oss.str();
  apf::writeVtkFiles(fname.c_str(), m);
  i++;
}

apf::Field* NestedAdjoint::run(double t_now, double t_old) {
  print_banner(t_now);
  solve(t_now, t_old);
  localize(t_now, t_old);
  auto bound = sum_contribs(u_error, p_error);
  auto e_nested = compute_error(u_error, p_error);
  print(" > |J(u)-J(uh)| < %.15e", bound);
  write_out();
  return nested_disc->set_error(e_nested);
}

NestedAdjoint* create_nested_adjoint(ParameterList const& p, Primal* pr) {
  return new NestedAdjoint(p, pr);
}

void destroy_nested_adjoint(NestedAdjoint* a) {
  delete a;
}

}
