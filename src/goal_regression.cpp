#include "goal_control.hpp"
#include "goal_functional.hpp"
#include "goal_regression.hpp"

namespace goal {

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<double>("tolerance", 0.0);
  p.set<double>("value", 0.0);
  return p;
}

void check_J_regression(ParameterList const& p, Functional* J) {
  if (! p.isSublist("regression")) return;
  auto rp = p.sublist("regression");
  rp.validateParameters(get_valid_params(), 0);
  auto tol = rp.get<double>("tolerance");
  auto val1 = rp.get<double>("value");
  auto val2 = J->get_value();
  auto err = std::abs(val2 - val1);
  GOAL_DEBUG_ASSERT(err < tol);
}

}
