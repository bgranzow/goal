#include <goal_control.hpp>
#include <goal_field.hpp>
#include <goal_state_fields.hpp>
#include <goal_traits.hpp>
#include <goal_workset.hpp>

#include "elast_ev_von_mises.hpp"

namespace elast {

template <typename EVALT, typename TRAITS>
VonMises<EVALT, TRAITS>::VonMises(RCP<goal::Field> u)
    : cauchy("Cauchy", u->get_tensor_ip_dl()),
      von_mises("Von Mises", u->get_scalar_ip_dl()) {
  /* populate the index dimensions */
  num_ips = u->get_num_elem_ips();
  num_dims = u->get_num_dims();

  /* populate the dependency structure for this evaluator. */
  this->addDependentField(cauchy);
  this->addEvaluatedField(von_mises);
  this->setName("Von Mises Stress");
}

template <typename EVALT, typename TRAITS>
void VonMises<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(cauchy, fm);
  this->utils.setFieldData(von_mises, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void VonMises<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  using Tensor = minitensor::Tensor<ScalarT>;
  Tensor sigma(3);
  sigma.fill(0.0);
  for (int elem = 0; elem < workset.size; ++elem) {
    for (int ip = 0; ip < num_ips; ++ip) {
      for (int i = 0; i < num_dims; ++i)
        for (int j = 0; j < num_dims; ++j)
          sigma(i, j) = cauchy(elem, ip, i, j);
      ScalarT s1 = (sigma(0,0)-sigma(1,1))*(sigma(0,0)-sigma(1,1));
      ScalarT s2 = (sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2));
      ScalarT s3 = (sigma(2,2)-sigma(0,0))*(sigma(2,2)-sigma(0,0));
      ScalarT s4 = sigma(0,1)*sigma(0,1);
      ScalarT s5 = sigma(1,2)*sigma(1,2);
      ScalarT s6 = sigma(2,0)*sigma(2,0);
      ScalarT s7 = 0.5*(s1+s2+s3+6.0*(s4+s5+s6));
      von_mises(elem, ip) = std::sqrt(s7);
    }
  }
}

template class VonMises<goal::Traits::Residual, goal::Traits>;
template class VonMises<goal::Traits::Jacobian, goal::Traits>;

}  /* namespace elast */
