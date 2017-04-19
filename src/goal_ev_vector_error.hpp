#ifndef GOAL_EV_VECTOR_ERROR_HPP
#define GOAL_EV_VECTOR_ERROR_HPP

/** \file goal_ev_vector_error.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
class Indexer;
class SolutionInfo;
/** \endcond */

/** \brief Fill in the global error field for a vector PDE.
  * \details This will scatter local contributions to the global error field
  * based on the intermediate residual evaluations. Element contributions
  * to the error field are considered for elements local to a workset.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * resid             | (Elem, Node)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * op                | (Dummy)
  */
template <typename EVALT, typename TRAITS>
class VectorError : public PHX::EvaluatorWithBaseImpl<TRAITS>,
                    public PHX::EvaluatorDerived<EVALT, TRAITS> {
 public:
  /** \cond */
  typedef typename TRAITS::SetupData SetupData;
  typedef typename TRAITS::PreEvalData PreEvalData;
  typedef typename TRAITS::PostEvalData PostEvalData;
  typedef typename TRAITS::EvalData EvalData;
  typedef typename EVALT::ScalarT ScalarT;
  /** \endcond */

  /** \brief Construct the evaluator.
    * \param u The primal solution vector DOF \ref goal::Field.
    * \param e The corresponding vector error \ref goal::Field.
    * \param i The relevant \ref goal::Indexer. */
  VectorError(RCP<Field> u, RCP<Field> e, RCP<Indexer> i);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Grab the linear algebra data. */
  void preEvaluate(PreEvalData info);

  /** \brief Fill in the local error field. */
  void evaluateFields(EvalData workset);
 
 private:
  RCP<Field> u;
  RCP<Field> e;
  RCP<Indexer> indexer;
  RCP<SolutionInfo> info;
  PHX::MDField<const ScalarT, Elem, Node> resid;
};

} /* namepsace goal */

#endif
