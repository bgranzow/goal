#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_linear_solve.hpp"

namespace goal {

typedef Tpetra::MultiVector<ST, LO, GO, KNode> MV;
typedef Tpetra::Operator<ST, LO, GO, KNode> OP;
typedef Tpetra::RowMatrix<ST, LO, GO, KNode> RM;
typedef Belos::LinearProblem<ST, MV, OP> LinearProblem;
typedef Belos::SolverManager<ST, MV, OP> Solver;
typedef Belos::BlockGmresSolMgr<ST, MV, OP> GmresSolver;
typedef Tpetra::Operator<ST, LO, GO, KNode> Prec;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<int>("krylov size", 0);
  p.set<int>("max iters", 0);
  p.set<double>("tolerance", 0.0);
  p.set<int>("output frequency", 0);
  p.set<int>("nonlinear max iters", 0);
  p.set<double>("nonlinear tolerance", 0.0);
  p.sublist("multigrid");
  return p;
}

static ParameterList get_belos_params(ParameterList const& in) {
  ParameterList p;
  int max_iters = in.get<int>("max iters");
  int krylov = in.get<int>("krylov size");
  double tol = in.get<double>("tolerance");
  p.set<int>("Block Size" , 1);
  p.set<int>("Num Blocks", krylov);
  p.set<int>("Maximum Iterations", max_iters);
  p.set<double>("Convergence Tolerance", tol);
  p.set<std::string>("Orthogonalization", "DGKS");
  if (in.isType<int>("output frequency")) {
    int f = in.get<int>("output frequency");
    p.set<int>("Verbosity", 33);
    p.set<int>("Output Style", 1);
    p.set<int>("Output Frequency", f);
  }
  return p;
}

static RCP<Solver> build_solver(
    ParameterList const& in,
    RCP<MatrixT> A,
    RCP<VectorT> x,
    RCP<VectorT> b,
    Disc* d) {
  Teuchos::ParameterList mg_params(in.sublist("multigrid"));
  auto belos_params = get_belos_params(in);
  auto AA = (RCP<OP>)A;
  auto coords = d->get_coords();
  auto P = MueLu::CreateTpetraPreconditioner(AA, mg_params, coords);
  auto problem = rcp(new LinearProblem(A, x, b));
  problem->setLeftPrec(P);
  problem->setProblem();
  return rcp(new GmresSolver(problem, rcpFromRef(belos_params)));
}

void solve(
    ParameterList const& in,
    RCP<MatrixT> A,
    RCP<VectorT> x,
    RCP<VectorT> b,
    Disc* d) {
  in.validateParameters(get_valid_params(), 0);
  auto solver = build_solver(in, A, x, b, d);
  auto dofs = solver->getProblem().getRHS()->getGlobalLength();
  print(" > linear system: num dofs %zu", dofs);
  auto t0 = time();
  solver->solve();
  auto t1 = time();
  auto iters = solver->getNumIters();
  print(" > linear system: solved in %d iterations", iters);
  if (iters >= in.get<int>("max iters"))
    print(" >  but solve was incomplete! continuing anyway...");
  print(" > linear system: solved in %f seconds", t1 - t0);
}

}
