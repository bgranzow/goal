#ifndef GOAL_EV_SCATTER_SCALAR_HPP
#define GOAL_EV_SCATTER_SCALAR_HPP

/** \file goal_ev_scatter_scalar.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_data_types.hpp"
#include "goal_dimension.hpp"
#include "goal_traits.hpp"

namespace goal {

/** \cond */
using Teuchos::RCP;

class Field;
class Indexer;
class SolutionInfo;

template <typename EvalT, typename Traits>
class ScatterScalar;
/** \endcond */


/** \brief Fill in the global residual vector from a scalar PDE.
  * \details This will scatter local contributions to the global residual
  * vector based on the intermediate physics evaluations. Element
  * contributions to the residual vector are considered for elements local
  * to a workset.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * resid             | (Elem, Node)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * op                | (Dummy)
  */
template <typename TRAITS>
class ScatterScalar<goal::Traits::Residual, TRAITS>
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
    * \param f The DOF \ref goal::Field corresponding to the scalar PDE.
    * \param i The relevant \ref goal::Indexer object.
    * \details The final unused boolean argument is to maintain interface
    * consistency with the Jacobian specialization of this class. */
  ScatterScalar(RCP<Field> f, RCP<Indexer> i, bool);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Grab the linear algebra data structures.
    * \param info The PreEvalData structure (\ref goal::SolutionInfo). */
  void preEvaluate(PreEvalData info);

  /** \brief Sum contributions to the residual vector. */
  void evaluateFields(EvalData workset);

 private:
  int num_nodes;
  RCP<Field> field;
  RCP<Indexer> indexer;
  RCP<SolutionInfo> info;
  PHX::MDField<const ScalarT, Elem, Node> resid;
};

/** \brief Fill in multiple global linear algebra data structures.
  * \details This will scatter local contributions to global linear algebra
  * data structures based on the intermediate physics evaluations and the
  * derivative information propogated through these evaluations with the
  * forward automatic differentiation data types. Local contributions to
  * both the Jacobian matrix and the residual vector are computed via this
  * evaluator. If the constructor argument `adj` is specified, then the
  * Jacobian matrix is filled in transpose form.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * resid             | (Elem, Node)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * op                | (Dummy)
  */
template <typename TRAITS>
class ScatterScalar<goal::Traits::Jacobian, TRAITS>
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
    * \param f The DOF \ref goal::Field corresponding to the scalar PDE.
    * \param i The relevant \ref goal::Indexer object.
    * \param adj True if the transpose of the Jacobian should be scattered. */
  ScatterScalar(RCP<Field> f, RCP<Indexer> i, bool adj);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Grab the linear algebra data structures.
    * \param info The PreEvalData structure (\ref goal::SolutionInfo). */
  void preEvaluate(PreEvalData info);

  /** \brief Sum contributions to the residual vector. */
  void evaluateFields(EvalData workset);

 private:
  void scatter_primal(
      RCP<Vector> R, RCP<Matrix> dRdu, int field_idx, EvalData workset);
  void scatter_adjoint(
      RCP<Vector> R, RCP<Matrix> dRduT, int field_idx, EvalData workset);

  int num_nodes;
  int num_dofs;
  int num_total_dofs;
  bool is_adjoint;
  RCP<Field> field;
  RCP<Indexer> indexer;
  RCP<SolutionInfo> info;
  PHX::MDField<const ScalarT, Elem, Node> resid;
};

}  // namespace goal

#endif
