#ifndef GOAL_EV_DIRICHLET_BCS_HPP
#define GOAL_EV_DIRICHLET_BCS_HPP

/** \file goal_ev_dirichlet_bcs.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"
#include "goal_traits.hpp"

namespace goal {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;

class Field;
class Indexer;
class SolutionInfo;

template <typename EVALT, typename TRAITS>
class DirichletBCs;
/** \endcond */

/** \brief Apply Dirichlet boundary conditions to the right hand side vector.
  * \details This evaluator will modify the right hand side of a linear system
  * to apply the appropriate boundary conditions for either the primal or
  * dual problem, as specified by this evaluator's constructor. For more
  * information about Dirichlet boundary conditions, please see the \ref DBCs
  * page. */
template <typename TRAITS>
class DirichletBCs<goal::Traits::Residual, TRAITS>
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
    * list should be a Tuechos::Array<std::string> of the form:
    * - [ dof field idx, dof field component, node set name, bc value ].
    * \param i The relevant \ref goal::Indexer for DOF information.
    * \param adj True if the dual boundary conditions should be applied. */
  DirichletBCs(RCP<const ParameterList> p, RCP<Indexer> i, bool adj);

  /** \brief Collect the linear algebra containers.
    * \param i The pre-evaluation data ( \ref goal::SolutionInfo ). */
  void preEvaluate(PreEvalData i);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
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
  * constructor. For more information about Dirichlet boundary conditions,
  * please see the \ref DBCs page. */
template <typename TRAITS>
class DirichletBCs<goal::Traits::Jacobian, TRAITS>
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
    * list should be a Tuechos::Array<std::string> of the form:
    * - [ dof field idx, dof field component, node set name, bc value ].
    * \param i The relevant \ref goal::Indexer for DOF information.
    * \param adj True if the dual boundary conditions should be applied. */
  DirichletBCs(RCP<const ParameterList> p, RCP<Indexer> i, bool adj);

  /** \brief Collect the linear algebra containers.
    * \param i The pre-evaluation data ( \ref goal::SolutionInfo ). */
  void preEvaluate(PreEvalData i);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  bool is_adjoint;
  RCP<Indexer> indexer;
  RCP<const ParameterList> params;
  RCP<SolutionInfo> info;
  void apply_bc(EvalData workset, Teuchos::Array<std::string> const& a);
};

}  // namespace goal

#endif
