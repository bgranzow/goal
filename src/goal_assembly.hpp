#include <Teuchos_RCP.hpp>

/** \file goal_assembly.hpp */

namespace goal {

using Teuchos::RCP;

/** \cond */
struct Workset;
class Physics;
class SolutionInfo;
class Discretization;
/** \endcond */

/** \brief Assemble the volumetric fields.
  * \param ws The \ref goal::Workset to operate on.
  * \param p The \ref goal::Physics to use for assembly.
  * \param i The \ref goal::Indexer used for linear algebra indexing.
  * \param d The \ref goal::Discretization to operate on. */
template <typename EvalT>
void assemble_volumetric(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);

/** \brief Assemble the Neumann fields.
  * \param ws The \ref goal::Workset to operate on.
  * \param p The \ref goal::Physics to use for assembly.
  * \param i The \ref goal::Indexer used for linear algebra indexing.
  * \param d The \ref goal::Discretization to operate on. */
template <typename EvalT>
void assemble_neumann(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);

/** \brief Assemble the Dirichlet fields.
  * \param ws The \ref goal::Workset to operate on.
  * \param p The \ref goal::Physics to use for assembly.
  * \param i The \ref goal::Indexer used for linear algebra indexing.
  * \param d The \ref goal::Discretization to operate on. */
template <typename EvalT>
void assemble_dirichlet(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);

/** \brief Compute the primal residual.
  * \param p The \ref goal::Physics to use for assembly.
  * \param i The \ref goal::Indexer used for linear algebra indexing.
  * \param d The \ref goal::Discretization to operate on.
  * \param t_current The time at the current step.
  * \param t_previous The time at the previous step. */
void compute_primal_residual(RCP<Physics> p, RCP<SolutionInfo> i,
    RCP<Discretization> d, const double t_current, const double t_previous);

/** \brief Compute the primal Jacobian.
  * \param p The \ref goal::Physics to use for assembly.
  * \param i The \ref goal::Indexer used for linear algebra indexing.
  * \param d The \ref goal::Discretization to operate on.
  * \param t_current The time at the current step.
  * \param t_previous The time at the previous step. */
void compute_primal_jacobian(RCP<Physics> p, RCP<SolutionInfo> i,
    RCP<Discretization> d, const double t_current, const double t_previous);

/** \brief Compute the dual Jacobian.
  * \param p The \ref goal::Physics to use for assembly.
  * \param i The \ref goal::Indexer used for linear algebra indexing.
  * \param d The \ref goal::Discretization to operate on.
  * \param t_current The time at the current step.
  * \param t_previous The time at the previous step. */
void compute_dual_jacobian(RCP<Physics> p, RCP<SolutionInfo> i,
    RCP<Discretization> d, const double t_current, const double t_previous);

/** \brief Compute the dual Jacobian.
  * \param p The \ref goal::Physics to use for assembly.
  * \param i The \ref goal::Indexer used for linear algebra indexing.
  * \param d The \ref goal::Discretization to operate on.
  * \param t_current The time at the current step.
  * \param t_previous The time at the previous step. */
void compute_error_residual(RCP<Physics> p, RCP<SolutionInfo> i,
    RCP<Discretization> d, const double t_current, const double t_previous);

}  /* namespace goal */
