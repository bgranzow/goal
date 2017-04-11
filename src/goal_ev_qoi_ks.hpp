#ifndef GOAL_EV_QOI_KS_HPP
#define GOAL_EV_QOI_KS_HPP

/** \file goal_ev_qoi_ks.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
/** \endcond */

/** \brief Compute the KS-functional of a scalar quantity \f$ q(u) \f$.
  * \details This evaluator computes the so-called KS-functional, which
  * attempts to approximate maximal values over the domain of a scalar
  * quantity \f$ g(u) \f$. The functional is given as:
  *
  * \f[
  * J(u) = \frac{1}{p} \ln \left[
  * \int_{\Omega} \exp{(p g(u))} \, \text{d} \Omega \right].
  * \f]
  *
  * and is implemented as:
  *
  * \f[
  * J(u) = m + \frac{1}{p} \ln \left[
  * \int_{\Omega} \exp{(p (g(u) - m))} \, \text{d} \Omega \right].
  * \f]
  *
  * for numerical reasons.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * g                 | (Elem, IP)
  * wdv               | (Elem, IP)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * J                 | (Elem)
  *
  * field descriptions:
  * - wdv, the differential volume (Jacobian) times the quadrature weight.
  * - g, the pre-computed scalar quantity evaluated at integrations points.
  * - J, the KS functional approximation to the maximal value \f$ \max g \f$.
  */
template <typename EVALT, typename TRAITS>
class QoIKS : public PHX::EvaluatorWithBaseImpl<TRAITS>,
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
    * \param u The relevant DOF \ref goal::Field \f$ u \f$.
    * \param n The name of the scalar quantity \f$ g(u) \f$.
    * \param p The scaling factor for the KS functional.
    * \param m The arbitrary numerical constant. */
  QoIKS(RCP<Field> u, std::string const& n, double p, double m);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

  /** \brief Apply global post-processing for the functional value. */
  void postEvaluate(PostEvalData data);

 private:
  int num_ips;
  double p;
  double m;
  double qoi_tmp;
  double qoi_value;
  PHX::MDField<double, Elem, IP> wdv;
  PHX::MDField<ScalarT, Elem, IP> g;
  PHX::MDField<ScalarT, Elem> J;
};

}  // namespace goal

#endif
