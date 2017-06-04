#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_control.hpp"
#include "goal_ev_qoi.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_sol_info.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;

template <typename TRAITS>
QoI<goal::Traits::Residual, TRAITS>::QoI(
    Indexer* i,
    Field* f,
    std::string const& n,
    int type) {

  auto dl = f->ent0_dl(type);
  qoi = PHX::MDField<const ScalarT, Ent>(n, dl);

  auto name = "QoI : " + n;
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));

  this->addDependentField(qoi);
  this->addEvaluatedField(op);
  this->setName(name);

  (void)i;
}

template <typename TRAITS>
void QoI<goal::Traits::Residual, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(qoi, fm);
  (void)d;
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
QoI<goal::Traits::Jacobian, TRAITS>::QoI(
    Indexer* i,
    Field* f,
    std::string const& n,
    int type)
    : indexer(i),
      info(0) {

  auto dl = f->ent0_dl(type);
  qoi = PHX::MDField<const ScalarT, Ent>(n, dl);
  num_dofs = indexer->get_num_total_dofs(type);

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
  auto dJdu = info->ghost->dJdu->get1dViewNonConst();
  std::vector<LO> rows(num_dofs);
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    indexer->get_ghost_lids(e, rows);
    for (int dof = 0; dof < num_dofs; ++dof) {
      LO row = rows[dof];
      auto val = qoi(elem).fastAccessDx(dof);
      dJdu[row] = val;
    }
  }
}

template class QoI<goal::Traits::Residual, goal::Traits>;
template class QoI<goal::Traits::Jacobian, goal::Traits>;

} // end namespace goal
