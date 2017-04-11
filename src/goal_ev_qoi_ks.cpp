#include "goal_ev_qoi_ks.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename EVALT, typename TRAITS>
QoIKS<EVALT, TRAITS>::QoIKS(
    RCP<Field> u, std::string const& name, double p_, double m_)
    : p(p_),
      m(m_),
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
  (void)workset;
}

template <typename EVALT, typename TRAITS>
void QoIKS<EVALT, TRAITS>::postEvaluate(PostEvalData info) {
  (void)info;
}

template class QoIKS<goal::Traits::Jacobian, goal::Traits>;

}  // namespace goal
