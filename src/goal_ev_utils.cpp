#include "goal_control.hpp"
#include "goal_ev_utils.hpp"
#include "goal_ev_gather_scalar.hpp"
#include "goal_ev_scalar_shape.hpp"
#include "goal_ev_interpolate_scalar.hpp"
#include "goal_ev_scatter_scalar.hpp"
#include "goal_ev_gather_vector.hpp"
#include "goal_ev_vector_shape.hpp"
#include "goal_ev_interpolate_vector.hpp"
#include "goal_ev_scatter_vector.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"

namespace goal {

using Teuchos::rcp;

template <typename EvalT>
void register_dof(RCP<Field> u, RCP<Indexer> i, FieldManager fm) {
  auto type = u->get_value_type();
  RCP<PHX::Evaluator<goal::Traits> > gather;
  RCP<PHX::Evaluator<goal::Traits> > shape;
  RCP<PHX::Evaluator<goal::Traits> > interpolate;
  if (type == SCALAR) {
    gather = rcp(new GatherScalar<EvalT, goal::Traits>(u, i));
    shape = rcp(new ScalarShape<EvalT, goal::Traits>(u));
    interpolate = rcp(new InterpolateScalar<EvalT, goal::Traits>(u));
  } else if (type == VECTOR) {
    gather = rcp(new GatherVector<EvalT, goal::Traits>(u, i));
    shape = rcp(new VectorShape<EvalT, goal::Traits>(u));
    interpolate = rcp(new InterpolateVector<EvalT, goal::Traits>(u));
  } else
    fail("cannot register dof");
  fm->registerEvaluator<EvalT>(gather);
  fm->registerEvaluator<EvalT>(shape);
  fm->registerEvaluator<EvalT>(interpolate);
}

template <typename EvalT>
void require_primal_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm) {
  auto type = u->get_value_type();
  RCP<PHX::Evaluator<goal::Traits> > scatter;
  if (type == SCALAR)
    scatter = rcp(new ScatterScalar<EvalT, goal::Traits>(u, i, false));
  else if (type == VECTOR)
    scatter = rcp(new ScatterVector<EvalT, goal::Traits>(u, i, false));
  else
    fail("cannot require primal scatter");
  fm->registerEvaluator<EvalT>(scatter);
  auto op = scatter->evaluatedFields()[0];
  fm->requireField<EvalT>(*op);
}

template <typename EvalT>
void require_adjoint_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm) {
  auto type = u->get_value_type();
  RCP<PHX::Evaluator<goal::Traits> > scatter;
  if (type == SCALAR)
    scatter = rcp(new ScatterScalar<EvalT, goal::Traits>(u, i, true));
  else if (type == VECTOR)
    scatter = rcp(new ScatterVector<EvalT, goal::Traits>(u, i, true));
  else
    fail("cannot require adjoint scatter");
  fm->registerEvaluator<EvalT>(scatter);
  auto op = scatter->evaluatedFields()[0];
  fm->requireField<EvalT>(*op);
}

template void register_dof<goal::Traits::Residual>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void register_dof<goal::Traits::Jacobian>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void require_primal_scatter<goal::Traits::Residual>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void require_primal_scatter<goal::Traits::Jacobian>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void require_adjoint_scatter<goal::Traits::Residual>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void require_adjoint_scatter<goal::Traits::Jacobian>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

}  // namespace goal
