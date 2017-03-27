#ifndef GOAL_EV_SCALAR_SHAPE_HPP
#define GOAL_EV_SCALAR_SHAPE_HPP

/** \file goal_ev_scalar_shape.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
/** \endcond */

/** \brief Compute the shape functions associated with a scalar field.
  * \details This evaluator fill in multidimensional arrays for the
  * shape functions and the shape function gradients for elements local
  * to the workset. This shape is fully determined by the \ref goal::Field
  * input, where \ref goal::Field::get_apf_basis and
  * \ref goal::Field::get_q_degree are used to fully determine these shape
  * functions.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * none              | N/A
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * wdv               | (Elem, IP)
  * shape             | (Elem, Node, IP)
  * grad_shape        | (Elem, Node, IP, Dim)
  *
  * field descriptions:
  * - wdv, The differential volume (Jacobian) of the element evaluated at
  * an integration point, weighted by the numerical quadrature weight at
  * that integration point.
  * - shape, The nodal shape functions associated with the input field
  * evaluated at integration points.
  * - grad_shape, The gradients of the nodal shape functions associated
  * with the input field evaluated at integration points. */
template <typename EVALT, typename TRAITS>
class ScalarShape : public PHX::EvaluatorWithBaseImpl<TRAITS>,
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
    * \param f The scalar \ref goal::Field to gather data from. */
  ScalarShape(RCP<Field> f);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_nodes;
  int num_ips;
  int num_dims;

  RCP<Field> field;

  PHX::MDField<double, Elem, IP> wdv;
  PHX::MDField<double, Elem, Node, IP> shape;
  PHX::MDField<double, Elem, Node, IP, Dim> grad_shape;
};

}  // namespace goal

#endif
