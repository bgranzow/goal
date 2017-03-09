#ifndef GOAL_EV_GATHER_SCALAR_HPP
#define GOAL_EV_GATHER_SCALAR_HPP

/** \file goal_ev_gather_scalar.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"
#include "goal_traits.hpp"

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
class Indexer;

template <typename EVALT, typename TRAITS>
class GatherScalar;
/** \endcond */

/** \brief Gather a scalar field locally from an existing APF field.
  * \details This will fill in a multidimensional array for the nodal
  * solution field that corresponds to elements local to the workset.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * none              | N/A
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * u                 | (Elem, Node)
 */
template <typename TRAITS>
class GatherScalar<goal::Traits::Residual, TRAITS>
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
    * \param f The scalar field to gather data from.
    * \param i The current relevant indexer for DOF information. */
  GatherScalar(RCP<Field> f, RCP<Indexer> i);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_nodes;
  RCP<Field> field;
  RCP<Indexer> indexer;
  PHX::MDField<ScalarT, Elem, Node> u;
};

/** \brief Gather a scalar field locally from an existing APF field.
  * \details This will fill in a multidimensional array for the nodal
  * solution field that corresponds to elements local to the workset.
  * This will also seed the elemental derivative array such that
  * \f$ \frac{\partial u^e_i}{\partial u^e_i} = \delta_{ij} \f$.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * none              | N/A
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * u                 | (Elem, Node)
  */
template <typename TRAITS>
class GatherScalar<goal::Traits::Jacobian, TRAITS>
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
    * \param f The scalar field to gather data from.
    * \param i The current relevant indexer for DOF information. */
  GatherScalar(RCP<Field> f, RCP<Indexer> i);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_nodes;
  RCP<Field> field;
  RCP<Indexer> indexer;
  PHX::MDField<ScalarT, Elem, Node> u;
};

}  // namespace goal

#endif
