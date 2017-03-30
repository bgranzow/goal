#include <goal_control.hpp>
#include <goal_field.hpp>
#include <Teuchos_ParameterList.hpp>

#include "elast_physics.hpp"

namespace elast {

using Teuchos::rcp;

static RCP<ParameterList> get_valid_params() {
  auto p = rcp(new ParameterList);
  p->sublist("dirichlet bcs");
  return p;
}

static void validate_params(RCP<const ParameterList> p) {
  assert(p->isType<std::string>("dirichlet bcs"));
  p->validateParameters(*get_valid_params(), 0);
}

Physics::Physics(
    RCP<const ParameterList> p, RCP<goal::Discretization> d)
    : goal::Physics(d) {
  params = p;
  validate_params(p);
}

Physics::~Physics() {
}

void Physics::set_primal() {
  is_primal = true;
  is_dual = false;
  is_error = false;
}

void Physics::set_dual() {
  is_primal = false;
  is_dual = true;
  is_error = false;
}

void Physics::set_error() {
  is_primal = false;
  is_dual = false;
  is_error = true;
}

}  // namespace elast
