#ifndef GOAL_EV_SCATTER_FUNCTIONAL
#define GOAL_EV_SCATTER_FUNCTIONAL

/** \file goal_ev_scatter_functional.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"
#include "goal_traits.hpp"

namespace goal {

/** \cond */
using Teuchos::RCP;

class Field;
class Indexer;

template <typename EVALT, typename TRAITS>
class ScatterFunctional;
/** \endcond */

/** \brief Fill in the functional derivative vector.
  * \details This evaluator will scatter local contributinos to the functional
  * derivative vector \f$ \frac{\partial J}{\partial u} \f$ based on element
  * level contributions to the function \f$ J(u) \f$ of the form:
  * \f[
  * J^e(u) = \int_{\Omega^e} j(u) \text{d} \Omega.
  * \f]
  * Notice that this class is only instantiated for the
  * \ref goal::Traits::Jacobian specialization.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * functional        | (Elem)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * op                | (Dummy)
  */
template <typename TRAITS>
class ScatterFunctional<goal::Traits::Jacobian, TRAITS>
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
    * \param f The relevant DOF \ref goal::Field.
    * \param i The relevant \ref goal::Indexer object.
    * \param qoi The name of the functional to be scattered. */
  ScatterFunctional(RCP<Field> f, RCP<Indexer> i, std::string const& qoi);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Grab the linear algebra data structures.
    * \param info The PreEvalData structure (\ref goal::SolutionInfo). */
  void preEvaluate(PreEvalData info);

  /** \brief Sum contributions to the functional derivative vector. */
  void evaluateFields(EvalData workset);

 private:
  int num_dofs;
  RCP<Field> field;
  RCP<Indexer> indexer;
  RCP<SolutionInfo> info;
  PHX::MDField<const ScalarT, Elem> functional;
};

}  /* namespace goal */

#endif
