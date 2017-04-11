#include <PCU.h>

#include "goal_ev_qoi_ks.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_solution_info.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename EVALT, typename TRAITS>
QoIKS<EVALT, TRAITS>::QoIKS(
    RCP<Field> u, std::string const& name, double p_, double m_)
    : p(p_),
      m(m_),
      qoi_tmp(0.0),
      qoi_val(0.0),
      g(name, u->get_scalar_ip_dl()),
      J("KS Functional", u->get_elem_scalar_dl()) {
  /* populate the index dimensions. */
  num_ips = u->get_num_elem_ips();

  /* set the dependency structure of this evaluator. */
  this->addDependentField(g);
  this->addEvaluatedField(J);
  this->setName("KS Functional");
}

template <typename EVALT, typename TRAITS>
void QoIKS<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(wdv, fm);
  this->utils.setFieldData(g, fm);
  this->utils.setFieldData(J, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void QoIKS<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  for (int elem = 0; elem < workset.size; ++elem) {
    J(elem) = 0.0;
    for (int ip = 0; ip < num_ips; ++ip) {
      qoi_tmp += std::exp(p*(g(elem, ip).val() - m)) * wdv(elem, ip);
      J(elem) += std::exp(p*(g(elem, ip) - m)) * wdv(elem, ip);
    }
  }
}

template <typename EVALT, typename TRAITS>
void QoIKS<EVALT, TRAITS>::postEvaluate(PostEvalData info) {
  auto dJdu = info.ghost->dJdu;
  PCU_Add_Doubles(&qoi_tmp, 1);
  double factor = 1.0/(p * qoi_tmp);
  dJdu->scale(factor);
  qoi_val = m + (1.0/p)*std::log(qoi_tmp);
}

template class QoIKS<goal::Traits::Jacobian, goal::Traits>;

}  // namespace goal
