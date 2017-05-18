#include "goal_control.hpp"
#include "goal_dbcs.hpp"
#include "goal_physics.hpp"
#include "goal_sol_info.hpp"

namespace goal {

void set_dbc_values(Physics* p, const double t) {
  (void)p;
  (void)t;
}

template <> void apply_primal_dbcs<goal::Traits::Residual>(
    Physics* p, SolInfo* i, const double t, bool c) {
  (void)p;
  (void)i;
  (void)t;
  (void)c;
}

template <> void apply_primal_dbcs<goal::Traits::Jacobian>(
    Physics* p, SolInfo* i, const double t, bool c) {
  (void)p;
  (void)i;
  (void)t;
  (void)c;
}

template <> void apply_dual_dbcs<goal::Traits::Jacobian>(
    Physics* p, SolInfo* i, const double t, bool c) {
  (void)p;
  (void)i;
  (void)t;
  (void)c;
}

} // end namespace goal
