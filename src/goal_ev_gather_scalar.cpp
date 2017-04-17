#include <apf.h>

#include "goal_control.hpp"
#include "goal_ev_gather_scalar.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename TRAITS>
GatherScalar<goal::Traits::Residual, TRAITS>::GatherScalar(
    RCP<Field> f, RCP<Indexer> i)
    : field(f), indexer(i), u(field->get_name(), field->get_dl()) {
  /* populate the index dimensions for this evaluator. */
  num_nodes = field->get_num_elem_nodes();
  assert(field->get_value_type() == SCALAR);

  /* populate the dependency structure of this evaluator. */
  this->addEvaluatedField(u);
  this->setName("Gather Scalar: " + field->get_name());
}

template <typename TRAITS>
void GatherScalar<goal::Traits::Residual, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(u, fm);
  (void)d;
}

template <typename TRAITS>
void GatherScalar<goal::Traits::Residual, TRAITS>::evaluateFields(
    EvalData workset) {
  auto f = field->get_apf_field();
  apf::NewArray<double> values;
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto fe = apf::createElement(f, e);
    apf::getScalarNodes(fe, values);
    for (int node = 0; node < num_nodes; ++node) u(elem, node) = values[node];
    apf::destroyElement(fe);
  }
}

template <typename TRAITS>
GatherScalar<goal::Traits::Jacobian, TRAITS>::GatherScalar(
    RCP<Field> f, RCP<Indexer> i)
    : field(f), indexer(i), u(field->get_name(), field->get_dl()) {
  /* populate the dimension information */
  num_nodes = field->get_num_elem_nodes();
  assert(field->get_value_type() == SCALAR);

  /* populate the dependency structure for this evaluator */
  this->addEvaluatedField(u);
  this->setName("Gather Scalar: " + field->get_name());
}

template <typename TRAITS>
void GatherScalar<goal::Traits::Jacobian, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(u, fm);
  (void)d;
}

template <typename TRAITS>
void GatherScalar<goal::Traits::Jacobian, TRAITS>::evaluateFields(
    EvalData workset) {
  auto f = field->get_apf_field();
  auto seed = field->get_seed_value();
  auto field_idx = field->get_associated_dof_idx();
  apf::NewArray<double> values;
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto fe = apf::createElement(f, e);
    apf::getScalarNodes(fe, values);
    for (int node = 0; node < num_nodes; ++node) {
      u(elem, node).val() = values[node];
      auto idx = indexer->get_elem_dof_offset(field_idx, node, 0);
      u(elem, node).fastAccessDx(idx) = seed;
    }
    apf::destroyElement(fe);
  }
}

template class GatherScalar<goal::Traits::Residual, goal::Traits>;
template class GatherScalar<goal::Traits::Jacobian, goal::Traits>;

}  /* namespace goal */
