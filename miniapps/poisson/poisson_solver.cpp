#include <apf.h>
#include <apfMesh.h>
#include <goal_assembly.hpp>
#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_ev_dirichlet_bcs.hpp>
#include <goal_field.hpp>
#include <goal_output.hpp>
#include <goal_solution_info.hpp>
#include <goal_linear_solvers.hpp>
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
  p->set<int>("num adapt cycles", 1);
  p->sublist("discretization");
  p->sublist("physics");
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
  void solve_primal();
  void solve_dual();
  void estimate_error();
  void adapt_mesh();
  RCP<const ParameterList> params;
  RCP<goal::Discretization> disc;
  RCP<poisson::Physics> physics;
  RCP<goal::SolutionInfo> info;
  RCP<goal::Output> output;
  int num_adapt_cycles;
};

Solver::Solver(RCP<const ParameterList> p) {
  params = p;
  validate_params(params);
  auto dp = rcpFromRef(params->sublist("discretization"));
  auto pp = rcpFromRef(params->sublist("physics"));
  auto op = rcpFromRef(params->sublist("output"));
  disc = rcp(new goal::Discretization(dp));
  physics = rcp(new poisson::Physics(pp, disc));
  output = rcp(new goal::Output(op, disc));
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
  physics->build_coarse_indexer();
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
  goal::add_to_fields(physics->get_u(), indexer, du);
  physics->destroy_model();
  physics->destroy_indexer();
}

void Solver::solve_dual() {
  goal::print("*** dual problem");
}

void Solver::estimate_error() {
  goal::print("*** error estimation");
}

void Solver::adapt_mesh() {
  goal::print("*** mesh adaptation");
}

void Solver::solve() {
  solve_primal();
  output->write(0);
}

} // namespace poisson


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
