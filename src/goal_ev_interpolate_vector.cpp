#include "goal_ev_interpolate_vector.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename EVALT, typename TRAITS>
InterpolateVector<EVALT, TRAITS>::InterpolateVector(RCP<Field> f)
    : field(f),
      shape(f->get_basis_name(), f->get_weight_dl()),
      grad_shape(f->get_grad_basis_name(), f->get_grad_weight_dl()),
      nodal(f->get_name(), f->get_dl()),
      u(f->get_name(), f->get_interpolated_dl()),
      grad_u(f->get_grad_name(), f->get_grad_interpolated_dl()) {

  /* populate the index dimensions for this evaluator. */
  num_nodes = field->get_num_elem_nodes();
  num_ips = field->get_num_elem_ips();
  num_dims = field->get_num_dims();
  assert(field->get_value_type() == VECTOR);

  /* populate the dependency structure of this evaluator. */
  this->addDependentField(shape);
  this->addDependentField(grad_shape);
  this->addDependentField(nodal);
  this->addEvaluatedField(u);
  this->addEvaluatedField(grad_u);
  this->setName("Interpolate: " + field->get_name());
}

template <typename EVALT, typename TRAITS>
void InterpolateVector<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(shape, fm);
  this->utils.setFieldData(grad_shape, fm);
  this->utils.setFieldData(nodal, fm);
  this->utils.setFieldData(u, fm);
  this->utils.setFieldData(grad_u, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void InterpolateVector<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  for (int elem = 0; elem < workset.size; ++elem) {
    for (int ip = 0; ip < num_ips; ++ip) {
      for (int i = 0; i < num_dims; ++i) {
        u(elem, ip, i) = nodal(elem, 0, i) * shape(elem, 0, ip, i);
        for (int node = 1; node < num_nodes; ++node)
          u(elem, ip, i) += nodal(elem, node, i) * shape(elem, node, ip, i);
        for (int j = 0; j < num_dims; ++j) {
          grad_u(elem, ip, i, j) =
              nodal(elem, 0, i) * grad_shape(elem, 0, ip, i, j);
          for (int node = 1; node < num_nodes; ++node)
            grad_u(elem, ip, i, j) +=
                nodal(elem, node, i) * grad_shape(elem, node, ip, i, j);
        }
      }
    }
  }
}

template class InterpolateVector<goal::Traits::Residual, goal::Traits>;
template class InterpolateVector<goal::Traits::Jacobian, goal::Traits>;

}  // namespace goal
