#ifndef GOAL_LINEAR_SOLVERS_HPP
#define GOAL_LINEAR_SOLVERS_HPP

/** \file goal_linear_solvers.hpp */

#include "goal_data_types.hpp"

namespace goal {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;
/** \endcond */

/** \brief Solve linear system with GMRES + ILU preconditioning.
  * \param p A parameter list that must contain:
  * - "linear max iters", int
  * - "linear krylov size", int
  * - "linear tolerance", double
  * \param A The matrix operator.
  * \param x The solution vector.
  * \param b The data vector. */
void solve_linear_system(
    RCP<const ParameterList> p, RCP<Matrix> A, RCP<Vector> x, RCP<Vector> b);

} // namespace goal

#endif
