#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>

#include "goal_control.hpp"
#include "goal_linear_solvers.hpp"

namespace goal {

using Teuchos::rcp;

typedef Tpetra::MultiVector<ST, LO, GO, KNode> MV;
typedef Tpetra::Operator<ST, LO, GO, KNode> OP;
typedef Tpetra::RowMatrix<ST, LO, GO, KNode> RM;
typedef Belos::LinearProblem<ST, MV, OP> LinearProblem;
typedef Belos::SolverManager<ST, MV, OP> Solver;
typedef Belos::BlockGmresSolMgr<ST, MV, OP> GmresSolver;
typedef Tpetra::Operator<ST, LO, GO, KNode> Prec;
typedef Ifpack2::Preconditioner<ST, LO, GO, KNode> IfpackPrec;


static RCP<ParameterList> get_valid_params() {
  auto p = rcp(new ParameterList);
  p->set<int>("linear max iters", 0);
  p->set<int>("linear krylov size", 0);
  p->set<double>("linear tolerance", 0.0);
  p->set<int>("linear output frequency", 0);
  return p;
}

static RCP<ParameterList> get_ifpack2_params() {
  auto p = rcp(new ParameterList);
  p->set("fact: drop tolerance", 0.0);
  p->set("fact: ilut level-of-fill", 1.0);
  return p;
}

static RCP<ParameterList> get_belos_params(RCP<const ParameterList> in) {
  auto p = rcp(new ParameterList);
  auto max_iters = in->get<int>("linear max iters");
  auto krylov = in->get<int>("linear krylov size");
  auto tol = in->get<double>("linear tolerance");
  p->set("Block Size", 1);
  p->set("Num Blocks", krylov);
  p->set("Maximum Iterations", max_iters);
  p->set("Convergence Tolerance", tol);
  p->set("Orthogonalization", "DGKS");
  if (in->isType<int>("linear output frequency")) {
    int f = in->get<int>("linear output frequency");
    p->set("Verbosity", 33);
    p->set("Output Style", 1);
    p->set("Output Frequency", f);
  }
  return p;
}

static RCP<Solver> build_ifpack2_solver(
    RCP<const ParameterList> in, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b) {
  auto ifpack2_params = get_ifpack2_params();
  auto belos_params = get_belos_params(in);
  Ifpack2::Factory factory;
  auto P = factory.create<RM>("ILUT", A);
  P->setParameters(*ifpack2_params);
  P->initialize();
  P->compute();
  auto problem = rcp(new LinearProblem(A, x, b));
  problem->setLeftPrec(P);
  problem->setProblem();
  auto solver = rcp(new GmresSolver(problem, belos_params));
  return solver;
}

static RCP<Solver> build_solver(
    RCP<const ParameterList> in, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b) {
  return build_ifpack2_solver(in, A, x, b);
}

void solve_linear_system(
    RCP<const ParameterList> in, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b) {
  auto t0 = time();
  in->validateParameters(*get_valid_params(), 0);
  print(" > linear system: num dofs %zu", x->getGlobalLength());
  auto solver = build_solver(in, A, x, b);
  solver->solve();
  auto iters = solver->getNumIters();
  auto t1 = time();
  if (iters >= in->get<int>("linear max iters"))
    print(
        " > linear solve failed to converge in %d iterations\n"
        " > continuing using the incomplete solve...",
        iters);
  else
    print(" > linear system solved in %d iterations", iters);
  print(" > linear system solved in %f seconds", t1 - t0);
}

}  // namespace goal
