#ifndef GOAL_EV_DUAL_SCALAR_WEIGHT_HPP
#define GOAL_EV_DUAL_SCALAR_WEIGHT_HPP

/** \file goal_ev_dual_scalar_weight.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
/** \endcond */

/** \brief Compute the dual weight times a partition of unity.
  * \details This evaluator fills in multidimensional arrays for the values
  * and gradients of a the dual weighting function \f$ (z^h - z^H) \f$
  * times a partition of \f$ \psi_i \f$, which is realized as linear Lagrange
  * shape functions. Here the evaluated weighting term at mesh vertices
  * is given as:
  *
  * \f[
  * w_i = (z^h - z^H) \psi_i,
  * \f]
  *
  * and the gradient of this term is given as:
  *
  * \f[
  * \nabla w_i = \nabla \left[ (z^h - z^H) \psi_i \right]
  * \f]
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * none              | N/A
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * w                 | (Elem, Node, IP)
  * grad_w            | (Elem, Node, IP, Dim)
  *
  * field descriptions:
  * - w, The dual weighting function times a partition of unity.
  * - grad_w, The gradient of w.
  * - note: Here the dimension tag Node is specialized to mesh vertices. */
template <typename EVALT, typename TRAITS>
class DualScalarWeight : public PHX::EvaluatorWithBaseImpl<TRAITS>,
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
    * \param z The dual \ref goal::Field on the coarse space \f$ V^H \f$.
    * \param z_fine The dual \ref goal::Field on the fine space \f$ V^h f$. */
  DualScalarWeight(RCP<Field> z, RCP<Field> z_fine);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  int num_vtx;
  int num_ips;
  int num_dims;
  RCP<Field> z;
  RCP<Field> z_fine;
  PHX::MDField<ScalarT, Elem, Node, IP> w;
  PHX::MDField<ScalarT, Elem, Node, IP, Dim> grad_w;
};

} /* namespace goal */

#endif
