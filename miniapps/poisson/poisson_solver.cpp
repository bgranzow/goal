#include <apf.h>
#include <apfMesh.h>
#include <ma.h>
#include <Teuchos_YamlParser_decl.hpp>
#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_output.hpp>
#include <goal_solution_info.hpp>


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
  RCP<goal::SolutionInfo> info;
  RCP<goal::Output> output;
  int num_adapt_cycles;
};

Solver::Solver(RCP<const ParameterList> p) {
  params = p;
  validate_params(params);
  auto dp = rcpFromRef(params->sublist("discretization"));
  disc = rcp(new goal::Discretization(dp));
}

void Solver::solve() {
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
