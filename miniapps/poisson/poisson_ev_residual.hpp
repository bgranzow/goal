#ifndef POISSON_EV_RESIDUAL_HPP
#define POISSON_EV_RESIDUAL_HPP

/** \file poisson_ev_residual.hpp */

#include <Phalanx_Evaluator_WithBaseImpl.hpp>
#include <Phalanx_Evaluator_Derived.hpp>
#include <goal_dimension.hpp>

/** \cond */
namespace goal {
class Field;
}
/** \endcond */

namespace poisson {

/** \cond */
using Teuchos::RCP;
/** \endcond */

/** \brief The residual evaluator for Poisson's equation.
  * \details This evaluator fills in multidimensional arrays for the residual
  * of Poisson's equation, corresponding to
  * \f[
  * \mathcal{R} (w, u) := ( \nabla w,  \nabla u) - (w, f).
  * \f]
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * wdv               | (Elem, IP)
  * w                 | (Elem, Node, IP)
  * grad_w            | (Elem, Node, IP, Dim)
  * grad_u            | (Elem, IP, Dim)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * resid             | (Elem, Node)
  *
  * field descriptions:
  * - wdv, The differential volume (Jacobian) of the element evaluated
  * at an integration point, multiplied by the numerical quadrature weight
  * at that integration point. Determined by \ref goal::ScalarShape.
  * - w, The nodal weighting function.
  * - grad_w, The gradient of the nodal weighting function.
  * - grad_u, The gradient of the solution function, as determined by
  * \ref goal::InterpolateScalar. */
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

  /** \brief Construct the evaluator for the forward or dual model.
    * \param u The scalar solution \ref goal::Field.
    * \param source The forcing function \f$ f \f$ for the problem. */
  Residual(RCP<goal::Field> u, std::string const& source);

  /** \brief Construct the evaluator for the error model.
    * \param u The scalar solution \ref goal::Field on the fine-scale space
    * \f$ V^h \f$.
    * \param z The scalar dual solution \ref goal::Field on the fine-scale
    * space \f$ V^h \f$.
    * \param source The forcing function \f$ f \f$ for the problem. */
  Residual(RCP<goal::Field> u, RCP<goal::Field> z, std::string const& source);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_nodes;
  int num_ips;
  int num_dims;
  RCP<goal::Field> sol;
  std::string ff;
  PHX::MDField<const double, Elem, IP> wdv;
  PHX::MDField<const double, Elem, Node, IP> w;
  PHX::MDField<const double, Elem, Node, IP, Dim> grad_w;
  PHX::MDField<const ScalarT, Elem, IP, Dim> grad_u;
  PHX::MDField<ScalarT, Elem, Node> resid;
};

}  /* namespace poisson */

#endif
