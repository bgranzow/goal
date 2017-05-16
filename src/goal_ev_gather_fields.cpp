#include <apf.h>

#include "goal_control.hpp"
#include "goal_ev_gather_fields.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename TRAITS>
GatherFields<goal::Traits::Residual, TRAITS>::GatherFields(
    Indexer* i, std::vector<Field*> const& f, int type)
    : indexer(i),
      fields(f) {
  num_fields = fields.size();
  num_nodes = fields[0]->get_num_nodes(type);
  u.resize(num_fields);
  for (int i = 0; i < num_fields; ++i) {
    u[i] = PHX::MDField<ScalarT, Ent, Node>(f[i]->name(), f[i]->dl(type));
    this->addEvaluatedField(u[i]);
  }
  this->setName("Gather Fields");
}

template <typename TRAITS>
void GatherFields<goal::Traits::Residual, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  for (int i = 0; i < num_fields; ++i)
    this->utils.setFieldData(u[i], fm);
  (void)d;
}

template <typename TRAITS>
void GatherFields<goal::Traits::Residual, TRAITS>::evaluateFields(
    EvalData workset) {
  apf::NewArray<double> values;
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    for (int field = 0; field < num_fields; ++field) {
      auto f = fields[field];
      auto apf_f = f->get_apf_field();
      auto fe = apf::createElement(apf_f, e);
      apf::getScalarNodes(fe, values);
      for (int n = 0; n < num_nodes; ++n)
        u[field](elem, n) = values[n];
      apf::destroyElement(fe);
    }
  }
}

template <typename TRAITS>
GatherFields<goal::Traits::Jacobian, TRAITS>::GatherFields(
    Indexer* i, std::vector<Field*> const& f, int type)
    : indexer(i),
      fields(f) {
  num_fields = fields.size();
  num_nodes = fields[0]->get_num_nodes(type);
  u.resize(num_fields);
  for (int i = 0; i < num_fields; ++i) {
    u[i] = PHX::MDField<ScalarT, Ent, Node>(f[i]->name(), f[i]->dl(type));
    this->addEvaluatedField(u[i]);
  }
  this->setName("Gather Fields");
}

template <typename TRAITS>
void GatherFields<goal::Traits::Jacobian, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  for (int i = 0; i < num_fields; ++i)
    this->utils.setFieldData(u[i], fm);
  (void)d;
}

template <typename TRAITS>
void GatherFields<goal::Traits::Jacobian, TRAITS>::evaluateFields(
    EvalData workset) {
  apf::NewArray<double> values;
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    for (int field = 0; field < num_fields; ++field) {
      auto f = fields[field];
      auto apf_f = f->get_apf_field();
      auto seed = f->get_seed_value();
      auto idx = f->get_associated_dof_idx();
      auto fe = apf::createElement(apf_f, e);
      apf::getScalarNodes(fe, values);
      for (int n = 0; n < num_nodes; ++n) {
        u[field](elem, n).val() = values[n];
        auto offset = indexer->get_elem_dof_offset(idx, n);
        u[field](elem, n).fastAccessDx(offset) = seed;
      }
      apf::destroyElement(fe);
    }
  }
}

template class GatherFields<goal::Traits::Residual, goal::Traits>;
template class GatherFields<goal::Traits::Jacobian, goal::Traits>;

} // end namespace goal
