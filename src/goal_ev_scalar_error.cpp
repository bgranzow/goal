#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_ev_scalar_error.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;

template <typename EVALT, typename TRAITS>
ScalarError<EVALT, TRAITS>::ScalarError(RCP<Field> u_, RCP<Field> e_)
    : u(u_),
      e(e_),
      resid(u->get_residual_name(), u->get_residual_PU_dl()) {
  /* make sure we're doing sane stuff. */
  assert(e->get_value_type() == SCALAR);
  assert(u->get_value_type() == SCALAR);

  /* populate the dependency structure for this evaluator. */
  auto name = "Error: " + u->get_name();
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));
  this->addDependentField(resid);
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename EVALT, typename TRAITS>
void ScalarError<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(resid, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void ScalarError<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  (void)workset;
}

template <typename EVALT, typename TRAITS>
void ScalarError<EVALT, TRAITS>::postEvaluate(PostEvalData info) {
  (void)info;
}

template class ScalarError<goal::Traits::Residual, goal::Traits>;

} /* namespace goal */
