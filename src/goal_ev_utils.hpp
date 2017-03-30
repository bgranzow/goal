#ifndef GOAL_EV_UTILS_HPP
#define GOAL_EV_UTILS_HPP

#include "goal_physics.hpp"

namespace goal {

template <typename EvalT>
void register_dof(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

template <typename EvalT>
void require_primal_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

template <typename EvalT>
void require_adjoint_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

}  // namespace goal

#endif
