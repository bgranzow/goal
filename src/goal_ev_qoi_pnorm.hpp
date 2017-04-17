#ifndef GOAL_EV_QOI_PNORM_HPP
#define GOAL_EV_QOI_PNORM_HPP

/** \file goal_ev_qoi_pnorm.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
/** \endcond */

/** \brief Compute the p-norm functional of a scalar quantity \f$ g(u) \f$.
  * \details This evaluator computes the so-called p-norm functional, which
  * attempts to approximate maximal values of the domain of a scalar
  * quantity \f$ g(u) \f$, and achieves the maximal value in the limit
  * \f$ p \to \infty \f$. The functional is given as:
  *
  * \f[
  * J(u) = \left[ \int_{\Omega} | g(u) |^p \, \text{d} \Omega
  * \right]^{\frac{1}{p}}
  * \f]
  *
  * and is implemented as:
  *
  * \f[
  * J(u) = m \left[ \int_{\Omega} | \frac{g(u)}{m} |^p \, \text{d} \Omega
  * \right]^{\frac{1}{p}}
  * \f]
  */
template <typename EVALT, typename TRAITS>
class QoIPNorm : public PHX::EvaluatorWithBaseImpl<TRAITS>,
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
    * \param p The scaling factor \f$ p \f$ for the p-norm functional.
    * \param m The arbitrary numerical constant \f$ m \f$. */
  QoIPNorm(RCP<Field> u, std::string const& n, double p, double m);

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
  double qoi_val;
  PHX::MDField<double, Elem, IP> wdv;
  PHX::MDField<ScalarT, Elem, IP> g;
  PHX::MDField<ScalarT, Elem> J;
};

}  /* namespace goal */

#endif
