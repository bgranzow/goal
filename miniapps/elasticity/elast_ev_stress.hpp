#ifndef ELAST_EV_STRESS_HPP
#define ELAST_EV_STRESS_HPP

/** \file elast_ev_stress.hpp */

#include <Phalanx_Evaluator_WithBaseImpl.hpp>
#include <Phalanx_Evaluator_Derived.hpp>
#include <goal_dimension.hpp>

/** \cond */
namespace goal {
class Field;
class StateFields;
}
/** \endcond */

namespace elast {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;
/** \endcond */

/** \brief Evaluates the Cauchy stress tensor.
  * \details This evaluator fills in a multidimensional array for the Cauchy
  * stress tensor for a linear isotropic constitutive model, given by:
  * \f[
  * \sigma = 2 \mu \epsilon + \lambda \text{tr}(\epsilon) I,
  * \f]
  * where \f$ \mu \f$ is the first Lame parameter given by:
  * \f[
  * \mu = \frac{E}{2(1 + \nu)},
  * \f]
  * \f$ \lambda \f$ is the second Lame parameter given by:
  * \f[
  * \lambda = \frac{E \nu}{(1 + \nu)(1 - 2 \nu)},
  * \f]
  * and \f$ \epsilon \f$ is the small strain tensor:
  * \f[
  * \epsilon = \frac12 (\nabla u + \nabla u^T).
  * \f]
  * Here \f$ E \f$ is the elastic modulus and \f$ \nu \f$ is Poisson's
  * ratio.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * grad_u            | (Elem, IP, Dim, Dim)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * cauchy            | (Elem, IP, Dim, Dim)
  *
  * field descriptions:
  * - grad_u, The displacement gradient evaluated at integration points.
  * - cauchy, The Cauchy stress tensor evaluated at integration points. */
template <typename EVALT, typename TRAITS>
class Stress : public PHX::EvaluatorWithBaseImpl<TRAITS>,
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
    * \param u The displacement \ref goal::Field.
    * \param s The \ref goal::StateFields structure where the Cauchy stress
    * tensor is stored.
    * \param mp The current element block's material parameters. Valid
    * parameters:
    * - "E", double, Elastic modulus.
    * - "nu", double, Poisson's ratio. */
  Stress(
      RCP<goal::Field> u,
      RCP<goal::StateFields> s,
      RCP<const ParameterList> mp);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  double E;
  double nu;
  int num_ips;
  int num_dims;
  RCP<goal::StateFields> states;
  PHX::MDField<const ScalarT, Elem, IP, Dim, Dim> grad_u;
  PHX::MDField<ScalarT, Elem, IP, Dim, Dim> cauchy;
};

}  /* namespace elast */

#endif
