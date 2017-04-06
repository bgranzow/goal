#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Ifpack2_Factory.hpp>

#include "goal_config.hpp"
#include "goal_control.hpp"
#include "goal_linear_solvers.hpp"

#ifdef Goal_MueLu
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#endif

namespace goal {

using Teuchos::rcp;

typedef Tpetra::MultiVector<ST, LO, GO, KNode> MV;
typedef Tpetra::Operator<ST, LO, GO, KNode> OP;
typedef Tpetra::RowMatrix<ST, LO, GO, KNode> RM;
typedef Belos::LinearProblem<ST, MV, OP> LinearProblem;
typedef Belos::SolverManager<ST, MV, OP> Solver;
typedef Belos::BlockCGSolMgr<ST, MV, OP> CGSolver;
typedef Belos::BlockGmresSolMgr<ST, MV, OP> GmresSolver;
typedef Tpetra::Operator<ST, LO, GO, KNode> Prec;
typedef Ifpack2::Preconditioner<ST, LO, GO, KNode> IfpackPrec;

static RCP<ParameterList> get_valid_params() {
  auto p = rcp(new ParameterList);
  p->set<int>("krylov size", 0);
  p->set<int>("maximum iterations", 0);
  p->set<double>("tolerance", 0.0);
  p->set<int>("output frequency", 0);
  p->sublist("multigrid");
  p->set<std::string>("method", "");
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
  int max_iters = in->get<int>("maximum iterations");
  int krylov = in->get<int>("krylov size");
  double tol = in->get<double>("tolerance");
  p->set<int>("Block Size" , 1);
  p->set<int>("Num Blocks", krylov);
  p->set<int>("Maximum Iterations", max_iters);
  p->set<double>("Convergence Tolerance", tol);
  p->set<std::string>("Orthogonalization", "DGKS");
  if (in->isType<int>("output frequency")) {
    int f = in->get<int>("output frequency");
    p->set<int>("Verbosity", 33);
    p->set<int>("Output Style", 1);
    p->set<int>("Output Frequency", f);
  }
  return p;
}

enum SolverType { CG, GMRES };

static RCP<Solver> build_ilu_solver(
    RCP<const ParameterList> in,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b,
    int type) {
  auto ilu_params = get_ifpack2_params();
  auto belos_params = get_belos_params(in);
  Ifpack2::Factory factory;
  auto P = factory.create<RM>("ILUT", A);
  P->setParameters(*ilu_params);
  P->initialize();
  P->compute();
  auto problem = rcp(new LinearProblem(A, x, b));
  problem->setLeftPrec(P);
  problem->setProblem();
  if (type == CG)
    return rcp(new CGSolver(problem, belos_params));
  else if (type == GMRES)
    return rcp(new GmresSolver(problem, belos_params));
  else
    return Teuchos::null;
}

static RCP<Solver> build_multigrid_solver(
    RCP<const ParameterList> in,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b,
    int type) {
  RCP<Solver> solver = Teuchos::null;
#ifdef Goal_MueLu
  Teuchos::ParameterList mg_params(in->sublist("multigrid"));
  auto belos_params = get_belos_params(in);
  auto AA = (RCP<OP>)A;
  auto P = MueLu::CreateTpetraPreconditioner(AA, mg_params);
  auto problem = rcp(new LinearProblem(A, x, b));
  problem->setLeftPrec(P);
  problem->setProblem();
  if (type == CG)
    solver = rcp(new CGSolver(problem, belos_params));
  else if (type == GMRES)
    solver = rcp(new GmresSolver(problem, belos_params));
#else
  (void)in;
  (void)A;
  (void)x;
  (void)b;
  (void)type;
  fail("calling multigrid preconditioner but Goal_MueLu=OFF!");
#endif
  return solver;
}

static void run_solver(
    RCP<const ParameterList> in,
    RCP<Solver> solver) {
  auto dofs = solver->getProblem().getRHS()->getGlobalLength();
  print(" > linear system: num dofs %zu", dofs);
  auto t0 = time();
  solver->solve();
  auto t1 = time();
  auto iters = solver->getNumIters();
  print(" > linear system: solved in %d iterations", iters);
  if (iters >= in->get<int>("maximum iterations"))
    print(" >  but solve was incomplete! continuing anyway...");
  print(" > linear system: solved in %f seconds", t1 - t0);
}

void solve_ilu_cg(
    RCP<const ParameterList> in,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b) {
  in->validateParameters(*get_valid_params(), 0);
  auto solver = build_ilu_solver(in, A, x, b, CG);
  run_solver(in, solver);
}

void solve_ilu_gmres(
    RCP<const ParameterList> in,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b) {
  in->validateParameters(*get_valid_params(), 0);
  auto solver = build_ilu_solver(in, A, x, b, GMRES);
  run_solver(in, solver);
}

void solve_multigrid_cg(
    RCP<const ParameterList> in,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b) {
  in->validateParameters(*get_valid_params(), 0);
  auto solver = build_multigrid_solver(in, A, x, b, CG);
  run_solver(in, solver);
}

void solve_multigrid_gmres(
    RCP<const ParameterList> in,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b) {
  in->validateParameters(*get_valid_params(), 0);
  auto solver = build_multigrid_solver(in, A, x, b, GMRES);
  run_solver(in, solver);
}

void solve_linear_system(
    RCP<const ParameterList> in,
    RCP<Matrix> A,
    RCP<Vector> x,
    RCP<Vector> b) {
  assert(in->isType<std::string>("method"));
  auto method = in->get<std::string>("method");
  assert( (method == "CG") || (method == "GMRES"));
  bool multigrid = false;
  if (in->isSublist("multigrid")) multigrid = true;
  if ((method == "CG") && (! multigrid))
    solve_ilu_cg(in, A, x, b);
  else if ((method == "GMRES") && (! multigrid))
    solve_ilu_gmres(in, A, x, b);
  else if ((method == "CG") && (multigrid))
    solve_multigrid_cg(in, A, x, b);
  else if ((method == "GMRES") && (multigrid))
    solve_multigrid_gmres(in, A, x, b);
  else
    fail("unkonwn linear solve combination");
}

}  // namespace goal
