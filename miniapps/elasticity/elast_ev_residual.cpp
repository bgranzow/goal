#include <goal_control.hpp>
#include <goal_field.hpp>
#include <goal_state_fields.hpp>
#include <goal_traits.hpp>
#include <goal_workset.hpp>

#include "elast_ev_residual.hpp"

namespace elast {

template <typename EVALT, typename TRAITS>
Residual<EVALT, TRAITS>::Residual(RCP<goal::Field> u)
    : wdv(u->get_wdv_name(), u->get_scalar_ip_dl()),
      grad_w(u->get_grad_basis_name(), u->get_grad_weight_dl()),
      cauchy("Cauchy", u->get_tensor_ip_dl()),
      resid(u->get_residual_name(), u->get_dl()) {
  /* populate index dimensions */
  num_nodes = u->get_num_elem_nodes();
  num_ips = u->get_num_elem_ips();
  num_dims = u->get_num_dims();

  /* populate the dependency structure for this evaluator */
  this->addDependentField(wdv);
  this->addDependentField(grad_w);
  this->addDependentField(cauchy);
  this->addEvaluatedField(resid);
  this->setName("Elastic Residual");
}

template <typename EVALT, typename TRAITS>
void Residual<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(wdv, fm);
  this->utils.setFieldData(grad_w, fm);
  this->utils.setFieldData(cauchy, fm);
  this->utils.setFieldData(resid, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void Residual<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  for (int elem = 0; elem < workset.size; ++elem) {
    for (int node = 0; node < num_nodes; ++node)
      for (int dim = 0; dim < num_dims; ++dim)
        resid(elem, node, dim) = ScalarT(0.0);
    for (int ip = 0; ip < num_ips; ++ip)
      for (int node = 0; node < num_nodes; ++node)
        for (int i = 0; i < num_dims; ++i)
          for (int j = 0; j < num_dims; ++j)
            resid(elem, node, i) += cauchy(elem, ip, i, j) *
                                    grad_w(elem, node, ip, i, j) *
                                    wdv(elem, ip);
  }
}

template class Residual<goal::Traits::Residual, goal::Traits>;
template class Residual<goal::Traits::Jacobian, goal::Traits>;

}  // namespace elast
