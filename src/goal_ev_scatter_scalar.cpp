#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_ev_scatter_scalar.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_solution_info.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;
using Teuchos::rcpFromRef;

template <typename TRAITS>
ScatterScalar<goal::Traits::Residual, TRAITS>::ScatterScalar(
    RCP<Field> f, RCP<Indexer> i, bool) {
  field = f;
  indexer = i;
  num_nodes = field->get_num_elem_nodes();
  auto name = "Scatter " + field->get_residual_name();
  resid = PHX::MDField<const ScalarT, Elem, Node>(
      field->get_residual_name(), field->get_dl());
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));
  this->addDependentField(resid);
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename TRAITS>
void ScatterScalar<goal::Traits::Residual, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(resid, fm);
  (void)d;
}

template <typename TRAITS>
void ScatterScalar<goal::Traits::Residual, TRAITS>::preEvaluate(PreEvalData i) {
  info = rcpFromRef(i);
  assert(info->ghost->R != Teuchos::null);
}

template <typename TRAITS>
void ScatterScalar<goal::Traits::Residual, TRAITS>::evaluateFields(
    EvalData workset) {
  auto R = info->ghost->R;
  auto field_idx = field->get_associated_dof_idx();
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    for (int node = 0; node < num_nodes; ++node) {
      LO row = indexer->get_ghost_lid(field_idx, e, node, 0);
      R->sumIntoLocalValue(row, resid(elem, node));
    }
  }
}

template <typename TRAITS>
ScatterScalar<goal::Traits::Jacobian, TRAITS>::ScatterScalar(
    RCP<Field> f, RCP<Indexer> i, bool adj) {
  field = f;
  indexer = i;
  is_adjoint = adj;
  num_nodes = field->get_num_elem_nodes();
  num_total_dofs = indexer->get_num_total_elem_dofs();
  auto name = "Scatter " + field->get_residual_name();
  resid = PHX::MDField<const ScalarT, Elem, Node>(
      field->get_residual_name(), field->get_dl());
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));
  this->addDependentField(resid);
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename TRAITS>
void ScatterScalar<goal::Traits::Jacobian, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(resid, fm);
  (void)d;
}

template <typename TRAITS>
void ScatterScalar<goal::Traits::Jacobian, TRAITS>::preEvaluate(PreEvalData i) {
  info = rcpFromRef(i);
  assert(info->ghost->R != Teuchos::null);
  assert(info->ghost->dRdu != Teuchos::null);
}

template <typename TRAITS>
void ScatterScalar<goal::Traits::Jacobian, TRAITS>::scatter_primal(
    RCP<Vector> R, RCP<Matrix> dRdu, int field_idx, EvalData workset) {
  using Teuchos::arrayView;
  std::vector<LO> cols(num_total_dofs);
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    indexer->get_ghost_lids(e, cols);
    auto c = arrayView(&cols[0], num_total_dofs);
    for (int node = 0; node < num_nodes; ++node) {
      auto v = resid(elem, node);
      auto view = arrayView(&(v.fastAccessDx(0)), num_total_dofs);
      LO row = indexer->get_ghost_lid(field_idx, e, node, 0);
      R->sumIntoLocalValue(row, v.val());
      dRdu->sumIntoLocalValues(row, c, view, num_total_dofs);
    }
  }
}

template <typename TRAITS>
void ScatterScalar<goal::Traits::Jacobian, TRAITS>::scatter_adjoint(
    RCP<Vector> R, RCP<Matrix> dRduT, int field_idx, EvalData workset) {
  using Teuchos::arrayView;
  std::vector<LO> cols(num_total_dofs);
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    indexer->get_ghost_lids(e, cols);
    for (int node = 0; node < num_nodes; ++node) {
      auto v = resid(elem, node);
      auto view = arrayView(&(v.fastAccessDx(0)), num_total_dofs);
      LO row = indexer->get_ghost_lid(field_idx, e, node, 0);
      R->sumIntoLocalValue(row, v.val());
      for (int dof = 0; dof < num_total_dofs; ++dof)
        dRduT->sumIntoLocalValues(
            cols[dof], arrayView(&row, 1), arrayView(&view[dof], 1));
    }
  }
}

template <typename TRAITS>
void ScatterScalar<goal::Traits::Jacobian, TRAITS>::evaluateFields(
    EvalData workset) {
  auto R = info->ghost->R;
  auto dRdu = info->ghost->dRdu;
  auto field_idx = field->get_associated_dof_idx();
  if (!is_adjoint)
    scatter_primal(R, dRdu, field_idx, workset);
  else
    scatter_adjoint(R, dRdu, field_idx, workset);
}

template class ScatterScalar<goal::Traits::Residual, goal::Traits>;
template class ScatterScalar<goal::Traits::Jacobian, goal::Traits>;

}  // namespace goal
