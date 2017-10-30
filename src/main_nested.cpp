#include <apf.h>
#include <Teuchos_YamlParameterListHelpers.hpp>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_functional.hpp"
#include "goal_mechanics.hpp"
#include "goal_nested_adjoint.hpp"
#include "goal_output.hpp"
#include "goal_primal.hpp"
#include "goal_regression.hpp"
#include "goal_states.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<int>("num steps", 0);
  p.set<double>("step size", 0.0);
  p.set<double>("initial time", 0.0);
  p.set<std::string>("adjoint mode", "");
  p.sublist("discretization");
  p.sublist("dirichlet bcs");
  p.sublist("traction bcs");
  p.sublist("mechanics");
  p.sublist("functional");
  p.sublist("primal linear algebra");
  p.sublist("adjoint linear algebra");
  p.sublist("output");
  p.sublist("regression");
  return p;
}

class Solver {
  public:
    Solver(const char* in);
    ~Solver();
    void solve();
  private:
    RCP<ParameterList> params;
    Disc* disc;
    Mechanics* mech;
    Primal* primal;
    Functional* functional;
    NestedAdjoint* adjoint;
    Output* output;
};

Solver::Solver(const char* in) {
  print("reading input: %s", in);
  params = rcp(new ParameterList);
  Teuchos::updateParametersFromYamlFile(in, params.ptr());
  params->validateParameters(get_valid_params(), 0);
  auto disc_params = params->sublist("discretization");
  auto mech_params = params->sublist("mechanics");
  auto out_params = params->sublist("output");
  disc = create_disc(disc_params);
  mech = create_mechanics(mech_params, disc);
  primal = create_primal(*params, mech);
  functional = create_functional(*params, primal);
  adjoint = create_nested_adjoint(*params, primal);
  output = create_output(out_params, disc);
}

Solver::~Solver() {
  destroy_output(output);
  destroy_functional(functional);
  destroy_nested_adjoint(adjoint);
  destroy_primal(primal);
  destroy_mechanics(mech);
  destroy_disc(disc);
}

void Solver::solve() {
  int num_steps = params->get<int>("num steps");
  double dt = params->get<double>("step size");
  double t_old = params->get<double>("initial time");
  double t_now = t_old + dt;
  for (int step = 1; step <= num_steps; ++step) {
    print("****** load step: %d", step);
    print("****** from time: %f", t_old);
    print("****** to time:   %f", t_now);
    disc->build_data();
    primal->build_data();
    primal->solve(t_now, t_old);
    functional->compute(t_now, t_old);
    functional->print_value();
    primal->destroy_data();
    disc->destroy_data();
    adjoint->build_data();
    auto e = adjoint->run(t_now, t_old);
    apf::destroyField(e);
    adjoint->destroy_data();
    output->write(t_now, 0);
    mech->get_states()->update();
    t_old = t_now;
    t_now += dt;
  }
  check_J_regression(*params, functional);
}

}

int main(int argc, char** argv) {
  goal::initialize();
  GOAL_DEBUG_ASSERT(argc == 2);
  const char* in = argv[1];
  { goal::Solver solver(in);
    solver.solve(); }
  goal::finalize();
}
