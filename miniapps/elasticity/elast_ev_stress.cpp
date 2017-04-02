#include <goal_control.hpp>
#include <goal_field.hpp>
#include <goal_state_fields.hpp>
#include <goal_traits.hpp>
#include <goal_workset.hpp>

#include "elast_ev_stress.hpp"

namespace elast {

static RCP<ParameterList> get_valid_params() {
  auto p = rcp(new ParameterList);
  p->set<double>("E", 0.0);
  p->set<double>("nu", 0.0);
  return p;
}

template <typename EVALT, typename TRAITS>
Stress<EVALT, TRAITS>::Stress(
    RCP<goal::Field> f, RCP<goal::StateFields> s, RCP<ParameterList> mp)
    : states(s),
      grad_u(f->get_grad_name(), f->get_grad_interpolated_dl()),
      cauchy("Cauchy", f->get_tensor_ip_dl()) {
  /* validate the material properties. */
  mp->validateParameters(*get_valid_params(), 0);

  /* populate material properies. */
  E = mp->get<double>("E");
  nu = mp->get<double>("nu");

  /* populate the dependency structure for this evaluator. */
  this->addDependentField(grad_u);
  this->addEvaluatedField(cauchy);
  this->setName("Elastic Stress");
}

template <typename EVALT, typename TRAITS>
void Stress<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(grad_u, fm);
  this->utils.setFieldData(cauchy, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void Stress<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  /* convenience */
  using Tensor = minitensor::Tensor<ScalarT>;

  /* declare tensor quantities */
  Tensor eps(num_dims), sigma(num_dims);
  Tensor I(minitensor::eye<ScalarT>(num_dims));

  /* define Lame material parameters */
  double mu = E / (2.0 * (1.0 + nu));
  double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

  /* loops for numerical integration */
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    for (int ip = 0; ip < num_ips; ++ip) {
      /* compute the small strain tensor */
      for (int i = 0; i < num_dims; ++i)
        for (int j = 0; j < num_dims; ++j)
          eps(i, j) = 0.5 * (grad_u(elem, ip, i, j) + grad_u(elem, ip, j, i));

      /* compute the Cauchy stress tensor */
      sigma = 2.0 * mu * eps + lambda * minitensor::trace(eps) * I;
      for (int i = 0; i < num_dims; ++i)
        for (int j = 0; j < num_dims; ++j) cauchy(elem, ip, i, j) = sigma(i, j);

      /* update the stress tensor state */
      states->set_tensor("cauchy", e, ip, sigma);
    }
  }
}

template class Stress<goal::Traits::Residual, goal::Traits>;
template class Stress<goal::Traits::Jacobian, goal::Traits>;

}  // namespace elast
