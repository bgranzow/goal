#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_ev_vector_error.hpp"
#include "goal_field.hpp"
#include "goal_solution_info.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;
using Teuchos::rcpFromRef;

template <typename EVALT, typename TRAITS>
VectorError<EVALT, TRAITS>::VectorError(
    RCP<Field> u_, RCP<Field> e_, RCP<Indexer> i)
    : u(u_),
      e(e_),
      indexer(i),
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
void VectorError<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(resid, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void VectorError<EVALT, TRAITS>::preEvaluate(PreEvalData i) {
  info = rcpFromRef(i);
  assert(info->owned->R != Teuchos::null);
}

template <typename EVALT, typename TRAITS>
void VectorError<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  (void)workset;
}

template class VectorError<goal::Traits::Residual, goal::Traits>;

} /* namespace goal */
