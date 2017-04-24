#ifndef GOAL_LINEAR_SOLVERS_HPP
#define GOAL_LINEAR_SOLVERS_HPP

/** \file goal_linear_solvers.hpp */

#include "goal_data_types.hpp"

namespace goal {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;

class Indexer;
/** \endcond */

/** \brief Conjugate Gradient solve with ILU preconditioning.
  * \param p A parameter-list with the following valid parameters:
  * - "krylov size", int, The size of the Krylov subspace.
  * - "maximum iterations", int, The maximum number of CG iterations.
  * - "tolerance", double, The CG convergence tolerance.
  * \param A The symmetric matrix operator.
  * \param x The linear system solution vector.
  * \param b The linear system data vector. */
void solve_ilu_cg(
    RCP<const ParameterList> p, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b);

/** \brief GMRES solve with ILU preconditioning.
  * \param p A parameter-list with the following valid parameters:
  * - "krylov size", int, The size of the Krylov subspace.
  * - "maximum iterations", int, The maximum number of GMRES iterations.
  * - "tolerance", double, The GMRES convergence tolerance.
  * \param A The matrix operator.
  * \param x The linear system solution vector.
  * \param b The linear system data vector. */
void solve_ilu_gmres(
    RCP<const ParameterList> p, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b);

/** \brief Conjugate Gradient solve with multigrid preconditioning.
  * \param p A parameterlist with the following valid parameters:
  * - "krylov size", int, The size of the Krylov subspace.
  * - "maximum iterations", int, The maximum number of CG iterations.
  * - "tolerance", double, The CG convergence tolerance.
  * - "multigrid", ParameterList, Valid MueLu parameters.
  * \param A The symmetric matrix operator.
  * \param x The linear system solution vector.
  * \param b The linear system data vector.
  * \param i Optionally pass the relevant \ref goal::Indexer to use the node
  * coordinates from the method \ref goal::Indexer::get_coords. This tends to
  * greatly improve preconditioner performance.
  * \details Trilinos must be compiled with MueLu to use this feature. */
void solve_multigrid_cg(
    RCP<const ParameterList> p, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b,
    RCP<Indexer> i = Teuchos::null);

/** \brief GMRES solve with multigrid preconditioning.
  * \param p A parameterlist with the following valid parameters:
  * - "krylov size", int, The size of the Krylov subspace.
  * - "maximum iterations", int, The maximum number of GMRES iterations.
  * - "tolerance", double, The GMRES convergence tolerance.
  * - "multigrid", ParameterList, Valid MueLu parameters.
  * \param A The symmetric matrix operator.
  * \param x The linear system solution vector.
  * \param b the linear system data vector.
  * \param i Optionally pass the relevant \ref goal::Indexer to use the node
  * coordinates from the method \ref goal::Indexer::get_coords. This tends to
  * greatly improve preconditioner performance.
  * \details Trilinos must be compiled with MueLu to use this feature. */
void solve_multigrid_gmres(
    RCP<const ParameterList> p, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b,
    RCP<Indexer> i = Teuchos::null);

/** \brief Solve a linear system iteratively.
  * \param p A parameterlist with the following valid parameters:
  * - "krylov size", int, The size of the Krylov subspace.
  * - "maximum iterations", int, The maximum number of iterations.
  * - "tolerance", double, The linear convergence tolerance.
  * - "method", std::string, The iterative method (CG, GMRES).
  * - "multigrid", ParameterList, Optional valid MueLu parameters.
  * \param A The matrix operator.
  * \param x The linear system solution vector.
  * \param b The linear system data vector.
  * \param i Optionally pass the relevant \ref goal::Indexer to use the node
  * coordinates from the method \ref goal::Indexer::get_coords. This tends to
  * greatly improve preconditioner performance. This is only used for
  * multigrid preconditioned solves.
  * \details If the parameterlist "multigrid" exists, then multigrid
  * preconditioning will be used. Otherwise, ILU preconditioning will
  * be defaulted. */
void solve_linear_system(
    RCP<const ParameterList> p, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b,
    RCP<Indexer> i = Teuchos::null);

} /* namespace goal */

#endif
