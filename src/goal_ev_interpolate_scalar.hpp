#ifndef GOAL_EV_INTERPOLATE_SCALAR_HPP
#define GOAL_EV_INTERPOLATE_SCALAR_HPP

/** \file goal_ev_interpolate_scalar.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
/** \endcond */

/** \brief Interpolate a nodal field to integration points.
  * \details This evaluator will fill in multidimensional arrays for the
  * values and gradients of nodal fields interpolated to the integration
  * points of an element.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * shape             | (Elem, Node, IP)
  * grad_shape        | (Elem, Node, IP, Dim)
  * nodal             | (Elem, Node)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * u                 | (Elem, IP)
  * grad_u            | (Elem, IP, Dim)
  *
  * field descriptions:
  * - shape, The input field's shape function values as determined by the
  * class \ref goal::ScalarShape.
  * - grad_shape, The input field's shape function gradient values as
  * determined by the class \ref goal::ScalarShape.
  * - nodal, The input field's nodal values as gathered by the class
  * goal::GatherScalar
  * - u, The values of the field interpolated to the integration points,
  * as determined by the method \ref goal::Field::get_q_degree.
  * - grad_u, The values of the gradient of the input field interpolated
  * to integration points, which are determined by the method
  * \ref goal::Field::get_q_degree. */
template <typename EVALT, typename TRAITS>
class InterpolateScalar
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
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
    * \param f The \ref goal::Field to interpolate. */
  InterpolateScalar(RCP<Field> f);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:

  int num_nodes;
  int num_ips;
  int num_dims;

  RCP<Field> field;

  PHX::MDField<const double, Elem, Node, IP> shape;
  PHX::MDField<const double, Elem, Node, IP, Dim> grad_shape;
  PHX::MDField<const ScalarT, Elem, Node> nodal;

  PHX::MDField<ScalarT, Elem, IP> u;
  PHX::MDField<ScalarT, Elem, IP, Dim> grad_u;
};

}  // namepsace goal

#endif
