#include <apf.h>

#include "goal_control.hpp"
#include "goal_ev_gather_vector.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename TRAITS>
GatherVector<goal::Traits::Residual, TRAITS>::GatherVector(
    RCP<Field> f, RCP<Indexer> i)
    : field(f),
      indexer(i),
      u(field->get_name(), field->get_dl()) {

  /* populate index dimensions for this evaluator. */
  assert(field->get_value_type() == VECTOR);
  num_nodes = field->get_num_elem_nodes();
  num_dims = field->get_num_dims();

  /* populate the dependecy structure of this evaluator. */
  this->addEvaluatedField(u);
  this->setName("Gather Vector: " + field->get_name());
}

template <typename TRAITS>
void GatherVector<goal::Traits::Residual, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(u, fm);
}

template <typename TRAITS>
void GatherVector<goal::Traits::Residual, TRAITS>::evaluateFields(
    EvalData workset) {
  auto f = field->get_apf_field();
  apf::NewArray<apf::Vector3> values;
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto fe = apf::createElement(f, e);
    apf::getVectorNodes(fe, values);
    for (int node = 0; node < num_nodes; ++node)
      for (int dim = 0; dim < num_dims; ++dim)
        u(elem, node, dim) = values[node][dim];
    apf::destroyElement(fe);
  }
}

template <typename TRAITS>
GatherVector<goal::Traits::Jacobian, TRAITS>::GatherVector(
    RCP<Field> f, RCP<Indexer> i)
    : field(f),
      indexer(i),
      u(field->get_name(), field->get_dl()) {

  /* populate index dimensions for this evaluator. */
  assert(field->get_value_type() == VECTOR);
  num_nodes = field->get_num_elem_nodes();
  num_dims = field->get_num_dims();

  /* populate the dependency structure for this evaluator. */
  this->addEvaluatedField(u);
  this->setName("Gather Vector: " + field->get_name());
}

template <typename TRAITS>
void GatherVector<goal::Traits::Jacobian, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(u, fm);
}

template <typename TRAITS>
void GatherVector<goal::Traits::Jacobian, TRAITS>::evaluateFields(
    EvalData workset) {
  auto f = field->get_apf_field();
  auto seed = field->get_seed_value();
  auto field_idx = field->get_associated_dof_idx();
  apf::NewArray<apf::Vector3> values;
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto fe = apf::createElement(f, e);
    apf::getVectorNodes(fe, values);
    for (int node = 0; node < num_nodes; ++node) {
      for (int dim = 0; dim < num_dims; ++dim) {
        u(elem, node, dim).val() = values[node][dim];
        auto idx = indexer->get_elem_dof_offset(field_idx, node, dim);
        u(elem, node, dim).fastAccessDx(idx) = seed;
      }
    }
    apf::destroyElement(fe);
  }
}

template class GatherVector<goal::Traits::Residual, goal::Traits>;
template class GatherVector<goal::Traits::Jacobian, goal::Traits>;

}  // namespace goal
