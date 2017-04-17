#ifndef ELAST_EV_RESIDUAL_HPP
#define ELAST_EV_RESIDUAL_HPP

/** \file elast_ev_residual.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>
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

/** \brief Evaluates the balance of linear momentum residual.
  * \details Assuming infinitesimal strains and given the Cauchy stress
  * tensor \f$ \sigma \f$, this evaluator fills in a multidimensional array
  * corresponding to the balance of linear momentum residual, given by the
  * the elemental integral:
  * \f[
  * R^e_i = \int_{\Omega^e} \sigma_{ij} w_{i,j} \; \text{d} \Omega,
  * \f]
  * where summation on \f$ i \f$ is not repeated.
  *
  * dependnet fields  | data layout
  * ----------------  | -----------
  * wdv               | (Elem, IP)
  * grad_w            | (Elem, Node, IP, Dim, Dim)
  * cauchy            | (Elem, IP, Dim, Dim)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * resid             | (Elem, Node, Dim)
  *
  * field descriptions:
  * - wdv, The elemental differential volume (Jacobian) at integration
  * points multiplied by the numerical integration weight.
  * - grad_w, The gradients of the nodal weighting functions evaluated
  * at integration points.
  * - cauchy, The Cauchy stress tensor evaluated at integration points.
  * - resid, Elemental contributions to the residual vector. */
template <typename EVALT, typename TRAITS>
class Residual : public PHX::EvaluatorWithBaseImpl<TRAITS>,
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
  Residual(RCP<goal::Field> u);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local mutlidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_nodes;
  int num_ips;
  int num_dims;
  PHX::MDField<const double, Elem, IP> wdv;
  PHX::MDField<const double, Elem, Node, IP, Dim, Dim> grad_w;
  PHX::MDField<const ScalarT, Elem, IP, Dim, Dim> cauchy;
  PHX::MDField<ScalarT, Elem, Node, Dim> resid;
};

}  /* namespace elast */

#endif
