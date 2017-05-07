#include <apf.h>
#include <apfMesh2.h>
#include <goal_assembly.hpp>
#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_ev_dirichlet_bcs.hpp>
#include <goal_error.hpp>
#include <goal_field.hpp>
#include <goal_indexer.hpp>
#include <goal_log.hpp>
#include <goal_output.hpp>
#include <goal_size_field.hpp>
#include <goal_solution_info.hpp>
#include <goal_linear_solvers.hpp>
#include <ma.h>
#include <Teuchos_YamlParameterListHelpers.hpp>

#include "elast_physics.hpp"

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Teuchos::ParameterList;

namespace elast {

static RCP<ParameterList> get_valid_params() {
  auto p = rcp(new ParameterList);
  p->sublist("discretization");
  p->sublist("physics");
  p->sublist("adaptation");
  p->sublist("linear algebra");
  p->sublist("output");
  return p;
}

static void validate_params(RCP<const ParameterList> p) {
  assert(p->isSublist("discretization"));
  assert(p->isSublist("physics"));
  assert(p->isSublist("linear algebra"));
  assert(p->isSublist("output"));
  p->validateParameters(*get_valid_params(), 0);
}

class Solver {
 public:
  Solver(RCP<const ParameterList> p);
  void solve();

 private:
  void ensure_quadratic();
  void solve_primal();
  void solve_dual();
  void estimate_error();
  void adapt_mesh();
  void solve_primal_only();
  void solve_adaptively();
  RCP<const ParameterList> params;
  RCP<goal::Discretization> disc;
  RCP<elast::Physics> physics;
  RCP<goal::SolutionInfo> info;
  RCP<goal::Output> output;
  RCP<goal::Log> log;
  bool should_adapt;
};

Solver::Solver(RCP<const ParameterList> p) {
  params = p;
  validate_params(params);
  auto dp = rcpFromRef(params->sublist("discretization"));
  auto pp = rcpFromRef(params->sublist("physics"));
  auto op = rcpFromRef(params->sublist("output"));
  disc = rcp(new goal::Discretization(dp));
  physics = rcp(new elast::Physics(pp, disc));
  output = rcp(new goal::Output(op, disc));
  log = rcp(new goal::Log(true, false));
  if (params->isSublist("adaptation")) should_adapt = true;
  else should_adapt = false;
}

static void zero_fields(RCP<Physics> p) {
  auto u = p->get_u();
  for (size_t i = 0; i < u.size(); ++i) {
    auto f = u[i]->get_apf_field();
    apf::zeroField(f);
  }
}

void Solver::solve_primal() {
  goal::print("*** primal problem");
  zero_fields(physics);
  physics->build_coarse_indexer(goal::STRIDED);
  physics->build_primal_model();
  goal::set_dbc_values(physics, 0.0);
  auto indexer = physics->get_indexer();
  info = rcp(new goal::SolutionInfo(indexer));
  goal::compute_primal_jacobian(physics, info, disc, 0.0, 0.0);
  auto dRdu = info->owned->dRdu;
  auto R = info->owned->R;
  auto du = info->owned->du;
  du->putScalar(0.0);
  R->scale(-1.0);
  auto lp = rcpFromRef(params->sublist("linear algebra"));
  goal::solve_linear_system(lp, dRdu, du, R, indexer);
  indexer->add_to_fields(physics->get_u(), du);
  goal::compute_primal_residual(physics, info, disc, 0.0, 0.0);
  log->pDOFs.push_back(du->getGlobalLength());
  physics->destroy_model();
  physics->destroy_indexer();
}

void Solver::solve_dual() {
  goal::print("*** dual problem");
  physics->build_enriched_data();
  physics->build_fine_indexer(goal::STRIDED);
  physics->build_dual_model();
  physics->enrich_state();
  auto indexer = physics->get_indexer();
  info = rcp(new goal::SolutionInfo(indexer, log));
  goal::compute_dual_jacobian(physics, info, disc, 0.0, 0.0);
  auto dRduT = info->owned->dRdu;
  auto dJdu = info->owned->dJdu;
  auto z = info->owned->z;
  z->putScalar(0.0);
  auto lp = rcpFromRef(params->sublist("linear algebra"));
  goal::solve_linear_system(lp, dRduT, z, dJdu, indexer);
  indexer->set_to_fields(physics->get_z_fine(), z);
  log->dDOFs.push_back(z->getGlobalLength());
  physics->destroy_model();
  physics->destroy_indexer();
}

static double dot_error(RCP<goal::SolutionInfo> i) {
  auto R = i->owned->R;
  auto z = i->owned->z;
  auto e = std::abs(R->dot(*z));
  goal::print(" > |J(u) - J(uH)| ~ %.15f", e);
  return e;
}

void Solver::estimate_error() {
  goal::print("*** error estimation");
  auto e = dot_error(info);
  physics->restrict_z_fine();
  physics->build_error_indexer(goal::STRIDED);
  physics->build_error_model();
  auto indexer = physics->get_indexer();
  info = rcp(new goal::SolutionInfo(indexer, log));
  goal::compute_error_residual(physics, info, disc, 0, 0);
  auto R = info->owned->R;
  indexer->set_to_fields(physics->get_e(), R);
  auto efields = physics->get_e();
  auto E_h = goal::sum_contributions(efields);
  auto B_h = goal::approx_upper_bound(efields);
  log->E_h.push_back(std::abs(E_h));
  log->B_h.push_back(B_h);
  std::cout << std::abs(e - std::abs(E_h)) << std::endl;
  physics->destroy_model();
  physics->destroy_indexer();
}

void Solver::adapt_mesh() {
  goal::print("*** mesh adaptation");
  static int scale = 1;
  auto m = disc->get_apf_mesh();
  auto efields = physics->get_e();
  auto i = goal::compute_indicators(efields);
  physics->destroy_enriched_data();
  physics->restrict_state();
  auto dofs = log->pDOFs.back();
  auto s = goal::get_iso_target_size(i, dofs * scale,  2);
  apf::destroyField(i);
  auto in = ma::configure(m, s);
  in->shouldRunPreParma = true;
  in->shouldRunPostParma = true;
  ma::adapt(in);
  apf::destroyField(s);
  scale *= 2;
}

void Solver::solve_primal_only() {
  solve_primal();
  output->write(0);
}

void Solver::solve_adaptively() {
  log->time.push_back(0.0);
  log->iter.push_back(0);
  solve_primal();
  solve_dual();
  estimate_error();
  adapt_mesh();
  output->write(0);
  log->print_summary();
}

void Solver::solve() {
  if (should_adapt) solve_adaptively();
  else solve_primal_only();
}

}  /* namespace elast */

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("Elasticity Solver!");
  assert(argc == 2);
  int status = EXIT_SUCCESS;
  try {
    const char* in = argv[1];
    goal::print("reading input file: %s", in);
    auto p = rcp(new ParameterList);
    Teuchos::updateParametersFromYamlFile(in, p.ptr());
    auto solver = rcp(new elast::Solver(p));
    solver->solve();
  } catch (std::exception const& ex) {
    goal::print("caught exception:");
    goal::print("%s", ex.what());
    status = EXIT_FAILURE;
  } catch (...) {
    goal::print("caught unknown exception");
    status = EXIT_FAILURE;
  }
  goal::finalize();
  return status;
}
