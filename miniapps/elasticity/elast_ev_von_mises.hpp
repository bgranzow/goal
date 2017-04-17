#ifndef ELAST_EV_VON_MISES_HPP
#define ELAST_EV_VON_MISES_HPP

/** \file elast_ev_von_mises.hpp */

#include <Phalanx_Evaluator_WithBaseImpl.hpp>
#include <Phalanx_Evaluator_Derived.hpp>
#include <goal_dimension.hpp>

/** \cond */
namespace goal {
class Field;
}
/** \endcond */

namespace elast {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;
/** \endcond */

/** \brief Evaluates the von Mises stress at integration points.
  * \details This evaluator fills in a multidimensional array for the
  * von Mises stress based on the Cauchy stress tensor, given by:
  * \f[
  * \sigma_{vm}^2 = \frac12 \left[
  * (\sigma_{11}-\sigma_{22})^2 +
  * (\sigma_{22}-\sigma_{33})^2 +
  * (\sigma_{33}-\sigma_{11})^2 +
  * 6 (\sigma_{12}^2 + \sigma_{23}^2 + \sigma_{31}^2)
  * \right]
  * \f]
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * cauchy            | (Elem, IP, Dim, Dim)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * von_mises         | (Elem, IP)
  *
  * field descriptions:
  * - cauchy, The Cauchy stress tensor evaluated at integration points.
  * - von_mises, The von Mises stress measure at integration points. */
template <typename EVALT, typename TRAITS>
class VonMises : public PHX::EvaluatorWithBaseImpl<TRAITS>,
                 public PHX::EvaluatorDerived<EVALT, TRAITS> {
 public:
  /** \cond */
  typedef typename TRAITS::SetupData SetupData;
  typedef typename TRAITS::PreEvalData PreEvalData;
  typedef typename TRAITS::PostEvalData PostEvalData;
  typedef typename TRAITS::EvalData EvalData;
  typedef typename EVALT::ScalarT ScalarT;
  using Elem = goal::Elem;
  using Node = goal::Node;
  using Dim = goal::Dim;
  using IP = goal::IP;
  /** \endcond */

  /** \brief Construct the evaluator.
    * \param u The displacement \ref goal::Field. */
  VonMises(RCP<goal::Field> u);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_ips;
  int num_dims;
  PHX::MDField<const ScalarT, Elem, IP, Dim, Dim> cauchy;
  PHX::MDField<ScalarT, Elem, IP> von_mises;
};

}  /* namespace elast */

#endif
