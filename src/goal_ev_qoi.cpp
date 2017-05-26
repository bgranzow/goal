#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_control.hpp"
#include "goal_ev_qoi.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_sol_info.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<int>("qoi index", 0);
  p.set<int>("ent type", 0);
  p.set<std::string>("name", "");
  p.set<goal::Field*>("field", 0);
  p.set<goal::Indexer*>("indexer", 0);
  return p;
}

template <typename TRAITS>
QoI<goal::Traits::Residual, TRAITS>::QoI(ParameterList const& p) {
  (void)p;
}

template <typename TRAITS>
void QoI<goal::Traits::Residual, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  (void)d;
  (void)fm;
}

template <typename TRAITS>
void QoI<goal::Traits::Residual, TRAITS>::preEvaluate(PreEvalData i) {
  (void)i;
}

template <typename TRAITS>
void QoI<goal::Traits::Residual, TRAITS>::evaluateFields(EvalData workset) {
  (void)workset;
}

template <typename TRAITS>
QoI<goal::Traits::Jacobian, TRAITS>::QoI(ParameterList const& p) {
  p.validateParameters(get_valid_params(), 0);

  auto type = p.get<int>("ent type");
  auto f = p.get<goal::Field*>("field");
  auto n = p.get<std::string>("name");
  auto dl = f->ent0_dl(type);

  qoi_idx = p.get<int>("qoi index");
  indexer = p.get<goal::Indexer*>("indexer");
  num_dofs = indexer->get_num_total_dofs(type);
  info = 0;

  qoi = PHX::MDField<const ScalarT, Ent>(n, dl);

  auto name = "QoI : " + n;
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));

  this->addDependentField(qoi);
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename TRAITS>
void QoI<goal::Traits::Jacobian, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(qoi, fm);
  (void)d;
}

template <typename TRAITS>
void QoI<goal::Traits::Jacobian, TRAITS>::preEvaluate(PreEvalData i) {
  info = i;
  GOAL_ALWAYS_ASSERT(Teuchos::nonnull(info->ghost->dJdu));
}

template <typename TRAITS>
void QoI<goal::Traits::Jacobian, TRAITS>::evaluateFields(EvalData workset) {
  auto dJdu = info->ghost->dJdu->getVector(qoi_idx);
  std::vector<LO> rows(num_dofs);
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    indexer->get_ghost_lids(e, rows);
    for (int dof = 0; dof < num_dofs; ++dof) {
      LO row = rows[dof];
      auto val = qoi(elem).fastAccessDx(dof);
      dJdu->sumIntoLocalValue(row, val);
    }
  }
}

template class QoI<goal::Traits::Residual, goal::Traits>;
template class QoI<goal::Traits::Jacobian, goal::Traits>;

} // end namespace goal
