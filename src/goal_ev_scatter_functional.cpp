#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_ev_scatter_functional.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_solution_info.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;
using Teuchos::rcpFromRef;

template <typename TRAITS>
ScatterFunctional<goal::Traits::Jacobian, TRAITS>::ScatterFunctional(
    RCP<Field> f, RCP<Indexer> i, std::string const& qoi) {
  field = f;
  indexer = i;
  num_dofs = indexer->get_num_total_elem_dofs();
  functional =
      PHX::MDField<const ScalarT, Elem>(qoi, field->get_elem_scalar_dl());
  auto name = "Scatter Functional";
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));
  this->addDependentField(functional);
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename TRAITS>
void ScatterFunctional<goal::Traits::Jacobian, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(functional, fm);
  (void)d;
}

template <typename TRAITS>
void ScatterFunctional<goal::Traits::Jacobian, TRAITS>::preEvaluate(
    PreEvalData i) {
  info = rcpFromRef(i);
  assert(info->ghost->dJdu != Teuchos::null);
}

template <typename TRAITS>
void ScatterFunctional<goal::Traits::Jacobian, TRAITS>::evaluateFields(
    EvalData workset) {
  auto dJdu = info->ghost->dJdu;
  std::vector<LO> rows;
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    indexer->get_ghost_lids(e, rows);
    for (int dof = 0; dof < num_dofs; ++dof) {
      LO row = rows[dof];
      auto val = functional(elem).fastAccessDx(dof);
      dJdu->sumIntoLocalValue(row, val);
    }
  }
}

template class ScatterFunctional<goal::Traits::Jacobian, goal::Traits>;

}  /* namespace goal */
