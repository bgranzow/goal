#include <PCU.h>
#include <apf.h>
#include <apfNumbering.h>

#include "goal_ev_qoi_scalar_point.hpp"
#include "goal_control.hpp"
#include "goal_field.hpp"
#include "goal_log.hpp"
#include "goal_traits.hpp"
#include "goal_indexer.hpp"
#include "goal_solution_info.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename EVALT, typename TRAITS>
QoIScalarPoint<EVALT, TRAITS>::QoIScalarPoint(
    RCP<Field> f, RCP<Indexer> i, std::string const& s)
    : field(f),
      indexer(i),
      set(s),
      J(0.0),
      u(f->get_name(), f->get_dl()),
      pw_u("Point-Wise " + f->get_name(), f->get_elem_scalar_dl()) {
  /* make sure we're doing sane stuff */
  assert(field->get_value_type() == SCALAR);

  /* grab the owned mesh vertex associated with the point x_0 */
  mesh = field->get_apf_mesh();
  auto idx = field->get_associated_dof_idx();
  auto nodes = indexer->get_node_set_nodes(set, idx);
  vtx = 0;
  if (nodes.size() > 0) {
    assert(nodes.size() == 1);
    vtx = nodes[0].entity;
  }

  /* set the dependency structure of this evaluator. */
  this->addDependentField(u);
  this->addEvaluatedField(pw_u);
  this->setName("Scalar Point Functional");
}

template <typename EVALT, typename TRAITS>
void QoIScalarPoint<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(u, fm);
  this->utils.setFieldData(pw_u, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void QoIScalarPoint<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  for (int elem = 0; elem < workset.size; ++elem)
    pw_u(elem) = 0.0;
}

template <typename EVALT, typename TRAITS>
void QoIScalarPoint<EVALT, TRAITS>::postEvaluate(PostEvalData info) {
  assert(Teuchos::nonnull(info.log));
  auto dJdu = info.ghost->dJdu->get1dViewNonConst();
  auto idx = field->get_associated_dof_idx();
  auto apf_field = field->get_apf_field();
  auto log = info.log;
  if (vtx && mesh->isOwned(vtx)) {
    LO row = indexer->get_ghost_lid(idx, vtx, 0, 0);
    dJdu[row] = 1.0;
    J = apf::getScalar(apf_field, vtx, 0);
  }
  PCU_Add_Doubles(&J, 1);
  print(" > J(u) ~ %.15f", J);
  log->Ju_h.push_back(J);
}

template class QoIScalarPoint<goal::Traits::Jacobian, goal::Traits>;

}  /* namespace goal */
