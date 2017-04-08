#ifndef GOAL_EV_CONDENSED_DBCS_HPP
#define GOAL_EV_CONDENSED_DBCS_HPP

/** \file goal_ev_condensed_dbcs.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"
#include "goal_traits.hpp"

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;
/** \endcond */

namespace goal {

/** \cond */
class Field;
class Indexer;
class SolutionInfo;

template <typename EVALT, typename TRAITS>
class CondensedDBCs;
/** \endcond */

/** \brief Apply Dirichlet boundary conditions to the residual vector.
  * \details This evaluator will modify the right hand side of a linear
  * system to apply the appropriate boundary conditions for either the 
  * primal or dual problem, as specified by this evaluator's constructor.
  * For more information about Dirichlet boundary conditions, please see
  * the \ref DBCs page. */
template <typename TRAITS>
class CondensedDBCs<goal::Traits::Residual, TRAITS>
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
      public PHX::EvaluatorDerived<goal::Traits::Residual, TRAITS> {
 public:
  /** \cond */
  typedef typename TRAITS::SetupData SetupData;
  typedef typename TRAITS::PreEvalData PreEvalData;
  typedef typename TRAITS::PostEvalData PostEvalData;
  typedef typename TRAITS::EvalData EvalData;
  typedef typename goal::Traits::Residual::ScalarT ScalarT;
  /** \endcond */

  /** \brief Construct the evaluator.
    * \param p The list of Dirichlet boundary conditions. Each entry in this
    * list should be a Teuchos::Array<std::string> of the form:
    * - [ dof field idx, dof field component, node set name, bc value ].
    * \param i The relevant \ref goal::Indexer for DOF information.
    * \param adj True if the dual boundary conditions should be applied. */
  CondensedDBCs(RCP<const ParameterList> p, RCP<Indexer> i, bool adj);

  /** \brief Collect the linear algebra containers.
    * \param i The pre-evaluation data (\ref goal::SolutionInfo). */
  void preEvaluate(PreEvalData i);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Modify the owned linear algebra data. */
  void evaluateFields(EvalData workset);

 private:
  bool is_adjoint;
  RCP<Indexer> indexer;
  RCP<const ParameterList> params;
  RCP<SolutionInfo> info;
  void apply_bc(EvalData workset, Teuchos::Array<std::string> const& a);
};

/** \brief Apply Dirichlet boundary conditions to the entire linear system.
  * \details This evaluator will modify both the left and right hand sides
  * of a linear system to apply the appropriate boundary conditions for
  * either the primal or dual problem, as specified by this evaluator's
  * constructor. In addition, the columns of the left hand side matrix are
  * condensed out. For more information about Dirichlet boundary conditions,
  * please see the \ref DBCs page */
template <typename TRAITS>
class CondensedDBCs<goal::Traits::Jacobian, TRAITS>
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
      public PHX::EvaluatorDerived<goal::Traits::Jacobian, TRAITS> {
 public:
  /** \cond */
  typedef typename TRAITS::SetupData SetupData;
  typedef typename TRAITS::PreEvalData PreEvalData;
  typedef typename TRAITS::PostEvalData PostEvalData;
  typedef typename TRAITS::EvalData EvalData;
  typedef typename goal::Traits::Jacobian::ScalarT ScalarT;
  /** \endcond */

  /** \brief Construct the evaluator.
    * \param p The list of Dirichlet boundary conditions. Each entry in this
    * list should be a Teuchos::Array<std::string> of the form:
    * - [ dof field idx, dof field component, node set name, bc value].
    * \param i The relevant \ref goal::Indexer for DOF information.
    * \param adj True if the dual boundary conditions should be applied. */
  CondensedDBCs(RCP<const ParameterList> p, RCP<Indexer> i, bool adj);

  /** \brief Collect the linear algebra containers.
    * \param i The pre-evaluation data (\ref goal::SolutionInfo). */
  void preEvaluate(PreEvalData i);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Modify the owned linear algebra data. */
  void evaluateFields(EvalData workset);

 private:
  bool is_adjoint;
  RCP<Indexer> indexer;
  RCP<const ParameterList> params;
  RCP<SolutionInfo> info;
};

}  // namespace goal

#endif
