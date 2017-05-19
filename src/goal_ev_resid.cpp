#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_control.hpp"
#include "goal_ev_resid.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_sol_info.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;

template <typename TRAITS>
Resid<goal::Traits::Residual, TRAITS>::Resid(
    Indexer* i, std::vector<Field*> const& f, int type, bool)
    : indexer(i),
      fields(f),
      info(0) {
  num_fields = fields.size();
  num_nodes = fields[0]->get_num_nodes(type);
  resids.resize(num_fields);
  for (int i = 0; i < num_fields; ++i) {
    auto n = fields[i]->resid_name();
    auto dl = fields[i]->dl(type);
    resids[i] = PHX::MDField<const ScalarT, Ent, Node>(n, dl);
    this->addDependentField(resids[i]);
  }
  PHX::Tag<ScalarT> op("Resid", rcp(new PHX::MDALayout<Dummy>(0)));
  this->addEvaluatedField(op);
  this->setName("Resid");
}

template <typename TRAITS>
void Resid<goal::Traits::Residual, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  for (int i = 0; i < num_fields; ++i)
    this->utils.setFieldData(resids[i], fm);
  (void)d;
}

template <typename TRAITS>
void Resid<goal::Traits::Residual, TRAITS>::preEvaluate(PreEvalData i) {
  info = i;
  GOAL_ALWAYS_ASSERT( Teuchos::nonnull(info->ghost->R) );
}

template <typename TRAITS>
void Resid<goal::Traits::Residual, TRAITS>::evaluateFields(
    EvalData workset) {
  auto R = info->ghost->R;
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    for (int node = 0; node < num_nodes; ++node) {
      for (int field = 0; field < num_fields; ++field) {
        LO row = indexer->get_ghost_lid(field, e, node);
        R->sumIntoLocalValue(row, resids[field](elem, node));
      }
    }
  }
}

template <typename TRAITS>
Resid<goal::Traits::Jacobian, TRAITS>::Resid(
    Indexer* i, std::vector<Field*> const& f, int type, bool is)
    : is_adjoint(is),
      indexer(i),
      fields(f),
      info(0) {
  num_fields = fields.size();
  num_nodes = fields[0]->get_num_nodes(type);
  num_dofs = indexer->get_num_total_dofs(type);
  resids.resize(num_fields);
  for (int i = 0; i < num_fields; ++i) {
    auto n = fields[i]->resid_name();
    auto dl = fields[i]->dl(type);
    resids[i] = PHX::MDField<const ScalarT, Ent, Node>(n, dl);
    this->addDependentField(resids[i]);
  }
  PHX::Tag<ScalarT> op("Resid", rcp(new PHX::MDALayout<Dummy>(0)));
  this->addEvaluatedField(op);
  this->setName("Resid");
}

template <typename TRAITS>
void Resid<goal::Traits::Jacobian, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  for (int i = 0; i < num_fields; ++i)
    this->utils.setFieldData(resids[i], fm);
  (void)d;
}

template <typename TRAITS>
void Resid<goal::Traits::Jacobian, TRAITS>::preEvaluate(PreEvalData i) {
  info = i;
  GOAL_ALWAYS_ASSERT( Teuchos::nonnull(info->ghost->R) );
  GOAL_ALWAYS_ASSERT( Teuchos::nonnull(info->ghost->dRdu) );
}

template <typename TRAITS>
void Resid<goal::Traits::Jacobian, TRAITS>::scatter_primal(
    RCP<Vector> R, RCP<Matrix> dRdu, EvalData workset) {
  using Teuchos::arrayView;
  std::vector<LO> cols(num_dofs);
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    indexer->get_ghost_lids(e, cols);
    auto c = arrayView(&cols[0], num_dofs);
    for (int node = 0; node < num_nodes; ++node) {
      for (int field = 0; field < num_fields; ++field) {
        auto v = resids[field](elem, node);
        auto view = arrayView(&(v.fastAccessDx(0)), num_dofs);
        LO row = indexer->get_ghost_lid(field, e, node);
        R->sumIntoLocalValue(row, v.val());
        dRdu->sumIntoLocalValues(row, c, view, num_dofs);
      }
    }
  }
}

template <typename TRAITS>
void Resid<goal::Traits::Jacobian, TRAITS>::scatter_adjoint(
    RCP<Vector> R, RCP<Matrix> dRduT, EvalData workset) {
  using Teuchos::arrayView;
  std::vector<LO> cols(num_dofs);
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    indexer->get_ghost_lids(e, cols);
    for (int node = 0; node < num_nodes; ++node) {
      for (int field = 0; field < num_fields; ++field) {
        auto v = resids[field](elem, node);
        auto view = arrayView(&(v.fastAccessDx(0)), num_dofs);
        LO row = indexer->get_ghost_lid(field, e, node);
        R->sumIntoLocalValue(row, v.val());
        for (int dof = 0; dof < num_dofs; ++dof)
          dRduT->sumIntoLocalValues(
              cols[dof], arrayView(&row, 1), arrayView(&view[dof], 1));
      }
    }
  }
}

template <typename TRAITS>
void Resid<goal::Traits::Jacobian, TRAITS>::evaluateFields(
    EvalData workset) {
  auto R = info->ghost->R;
  auto dRdu = info->ghost->dRdu;
  if (! is_adjoint) scatter_primal(R, dRdu, workset);
  else scatter_adjoint(R, dRdu, workset);
}

template class Resid<goal::Traits::Residual, goal::Traits>;
template class Resid<goal::Traits::Jacobian, goal::Traits>;

}
