#ifndef GOAL_EV_GATHER_VECTOR_HPP
#define GOAL_EV_GATHER_VECTOR_HPP

/** \file goal_ev_gather_vector.hpp */

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
class GatherVector;
/** \endcond */

/** \brief Gather a vector field locally from an existing APF field.
  * \details This will fill in a multidimensional array for the nodal
  * solution field that corresponds to elements local to the workset.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * none              | N/A
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * u                 | (Elem, Node, Dim)
  */
template <typename TRAITS>
class GatherVector<goal::Traits::Residual, TRAITS>
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
    * \param f The scalar \ref goal::Field to gather data from.
    * \param i The relevant \ref goal::Indexer for DOF information. */
  GatherVector(RCP<Field> f, RCP<Indexer> i);

  /** \brief Finalize the filed manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_nodes;
  int num_dims;
  RCP<Field> field;
  RCP<Indexer> indexer;
  PHX::MDField<ScalarT, Elem, Node, Dim> u;
};

/** \brief Gather a vector field locally from an existing APF field.
  * \details This will fill in a multidimensional array for the nodal
  * solution field that corresponds to elements local to the workst.
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
class GatherVector<goal::Traits::Jacobian, TRAITS>
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
    * \param f The vector \ref goal::Field to gather data from.
    * \param i The relevant \ref goal::Indexer for DOF information. */
  GatherVector(RCP<Field> f, RCP<Indexer> i);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_nodes;
  int num_dims;
  RCP<Field> field;
  RCP<Indexer> indexer;
  PHX::MDField<ScalarT, Elem, Node, Dim> u;
};

}  // namespace goal

#endif
