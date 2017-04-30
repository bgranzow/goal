#include <apf.h>
#include <apfMesh.h>
#include <goal_assembly.hpp>
#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_error.hpp>
#include <goal_ev_dirichlet_bcs.hpp>
#include <goal_log.hpp>
#include <goal_field.hpp>
#include <goal_indexer.hpp>
#include <goal_output.hpp>
#include <goal_solution_info.hpp>
#include <goal_linear_solvers.hpp>
#include <goal_size_field.hpp>
#include <ma.h>
#include <Teuchos_YamlParameterListHelpers.hpp>

#include "poisson_physics.hpp"

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Teuchos::ParameterList;

namespace poisson {

static RCP<ParameterList> get_valid_params() {
  auto p = rcp(new ParameterList);
  p->sublist("discretization");
  p->sublist("physics");
  p->sublist("linear algebra");
  p->sublist("adaptation");
  p->sublist("output");
  return p;
}

static RCP<ParameterList> get_valid_adapt_params() {
  auto p = rcp(new ParameterList);
  p->set<int>("num cycles", 0);
  p->set<int>("initial target elems", 0);
  p->set<double>("J exact", 0.0);
  return p;
}

static void validate_params(RCP<const ParameterList> p) {
  assert(p->isSublist("discretization"));
  assert(p->isSublist("physics"));
  assert(p->isSublist("linear algebra"));
  assert(p->isSublist("adaptation"));
  assert(p->isSublist("output"));
  auto ap = rcpFromRef(p->sublist("adaptation"));
  assert(ap->isType<int>("num cycles"));
  assert(ap->isType<int>("initial target elems"));
  assert(ap->isType<double>("J exact"));
  p->validateParameters(*get_valid_params(), 0);
  ap->validateParameters(*get_valid_adapt_params(), 0);
}

static RCP<goal::Log> create_log(RCP<const ParameterList> p) {
  auto J = p->get<double>("J exact");
  return rcp(new goal::Log(true, true, J));
}

class Solver {
 public:
  Solver(RCP<const ParameterList> p);
  void solve();

 private:
  void solve_primal();
  void solve_dual();
  void estimate_error();
  void adapt_mesh();
  RCP<const ParameterList> params;
  RCP<goal::Discretization> disc;
  RCP<poisson::Physics> physics;
  RCP<goal::SolutionInfo> info;
  RCP<goal::Output> output;
  RCP<goal::Log> log;
  int num_adapt_cycles;
  int target_init;
};

Solver::Solver(RCP<const ParameterList> p) {
  params = p;
  validate_params(params);
  auto dp = rcpFromRef(params->sublist("discretization"));
  auto pp = rcpFromRef(params->sublist("physics"));
  auto op = rcpFromRef(params->sublist("output"));
  auto ap = rcpFromRef(params->sublist("adaptation"));
  disc = rcp(new goal::Discretization(dp));
  physics = rcp(new poisson::Physics(pp, disc));
  output = rcp(new goal::Output(op, disc));
  log = create_log(ap);
  num_adapt_cycles = ap->get<int>("num cycles");
  target_init = ap->get<int>("initial target elems");
}

static void zero_fields(std::vector<RCP<goal::Field> > u) {
  for (size_t i = 0; i < u.size(); ++i) {
    auto f = u[i]->get_apf_field();
    apf::zeroField(f);
  }
}

void Solver::solve_primal() {
  goal::print("*** primal problem");
  auto u_fields = physics->get_u();
  zero_fields(u_fields);
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
  goal::solve_linear_system(lp, dRdu, du, R);
  indexer->add_to_fields(physics->get_u(), du);
  log->pDOFs.push_back(du->getGlobalLength());
  physics->destroy_model();
  physics->destroy_indexer();
}

void Solver::solve_dual() {
  goal::print("*** dual problem");
  physics->build_enriched_data();
  physics->build_fine_indexer(goal::STRIDED);
  physics->build_dual_model();
  auto indexer = physics->get_indexer();
  info = rcp(new goal::SolutionInfo(indexer, log));
  goal::compute_dual_jacobian(physics, info, disc, 0, 0);
  auto dRduT = info->owned->dRdu;
  auto dJdu = info->owned->dJdu;
  auto z = info->owned->z;
  auto lp = rcpFromRef(params->sublist("linear algebra"));
  goal::solve_linear_system(lp, dRduT, z, dJdu);
  indexer->set_to_fields(physics->get_z_fine(), z);
  log->dDOFs.push_back(z->getGlobalLength());
  physics->destroy_model();
  physics->destroy_indexer();
}

void Solver::estimate_error() {
  goal::print("*** error estimation");
  auto R = info->owned->R;
  auto z = info->owned->z;
  auto e = std::abs(R->dot(*z));
  goal::print(" > |J(u) - J(uH)| ~ %.15f", e);
  physics->restrict_z_fine();
  physics->build_error_indexer(goal::STRIDED);
  physics->build_error_model();
  goal::compute_error_residual(physics, info, disc, 0, 0);
  auto indexer = physics->get_indexer();
  indexer->set_to_fields(physics->get_e(), R);
  auto efields = physics->get_e();
  auto E_h = goal::sum_contributions(efields);
  auto B_h = goal::approx_upper_bound(efields);
  log->E_h.push_back(std::abs(E_h));
  log->B_h.push_back(B_h);
  assert(std::abs(e - std::abs(E_h)) < 1.0e-8);
  physics->destroy_model();
  physics->destroy_indexer();
}

void Solver::adapt_mesh() {
  goal::print("*** mesh adaptation");
  static int scale = 1;
  auto m = disc->get_apf_mesh();
  auto e = physics->get_e()[0]->get_apf_field();
  auto s = goal::get_iso_target_size(e, scale * target_init, 1);
  physics->destroy_enriched_data();
  auto in = ma::configure(m, s);
  in->shouldRunPostParma = true;
  ma::adapt(in);
  apf::destroyField(s);
  disc->update();
  scale *= 2;
}

void Solver::solve() {
  for (int i = 0; i <= num_adapt_cycles; ++i) {
    log->time.push_back(0.0);
    log->iter.push_back(i);
    solve_primal();
    solve_dual();
    estimate_error();
    output->write(i);
    adapt_mesh();
  }
  log->print_summary();
}

} /* namespace poisson */


int main(int argc, char** argv) {
  goal::initialize();
  goal::print("Poisson Solver!");
  assert(argc == 2);
  int status = EXIT_SUCCESS;
  try {
    const char* in = argv[1];
    goal::print("reading input file: %s", in);
    auto p = rcp(new ParameterList);
    Teuchos::updateParametersFromYamlFile(in, p.ptr());
    auto solver = rcp(new poisson::Solver(p));
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
