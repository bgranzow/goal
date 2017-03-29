#include <goal_control.hpp>
#include <goal_field.hpp>
#include <Teuchos_ParameterList.hpp>

#include "poisson_physics.hpp"

namespace poisson {

using Teuchos::rcp;

static RCP<ParameterList> get_valid_params() {
  auto p = rcp(new ParameterList);
  p->set<std::string>("forcing function", "");
  p->set<std::string>("point set", "");
  p->set<std::string>("functional type", "");
  p->sublist("dirichlet bcs");
  return p;
}

static void validate_params(RCP<const ParameterList> p) {
  assert(p->isType<std::string>("forcing function"));
  assert(p->isType<std::string>("functional type"));
  assert(p->isSublist("dirichlet bcs"));
  auto type = p->get<std::string>("functional type");
  if (type == "point-wise") assert(p->isType<std::string>("point set"));
  p->validateParameters(*get_valid_params(), 0);
}

Physics::Physics(
    RCP<const ParameterList> p, RCP<goal::Discretization> d)
    : goal::Physics(d) {
  params = p;
  validate_params(p);
  ff = params->get<std::string>("forcing function");
  functional_type = p->get<std::string>("functional type");
  if (functional_type == "point-wise")
    set = params->get<std::string>("point set");
  goal::FieldInfo info =
    {disc, "u", 1, 1, goal::SCALAR, goal::HIERARCHIC};
  u.push_back(rcp(new goal::Field(&info)));
  u[0]->set_associated_dof_idx(0);
  u[0]->set_dof_status(true);
  u[0]->set_seed(1.0);
  is_primal = false;
  is_dual = false;
  is_error = false;
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

void Physics::build_primal_volumetric(FieldManager fm) {
  set_primal();
  register_volumetric<Residual>(fm);
  register_volumetric<Jacobian>(fm);
}

void Physics::build_primal_dirichlet(FieldManager fm) {
  set_primal();
  register_dirichlet<Residual>(fm);
  register_dirichlet<Jacobian>(fm);
}

void Physics::build_dual_volumetric(FieldManager fm) {
  set_dual();
  register_volumetric<Jacobian>(fm);
}

void Physics::build_dual_dirichlet(FieldManager fm) {
  set_dual();
  register_dirichlet<Jacobian>(fm);
}

void Physics::build_error_volumetric(FieldManager fm) {
  set_error();
  register_volumetric<Residual>(fm);
}

}  // namespace poisson

template <typename EvalT>
void poisson::Physics::register_volumetric(goal::FieldManager fm) {
  (void)fm;
}

template <typename EvalT>
void poisson::Physics::register_dirichlet(goal::FieldManager fm) {
  (void)fm;
}

template void poisson::Physics::register_volumetric<goal::Traits::Residual>(goal::FieldManager fm);
template void poisson::Physics::register_volumetric<goal::Traits::Jacobian>(goal::FieldManager fm);
template void poisson::Physics::register_dirichlet<goal::Traits::Residual>(goal::FieldManager fm);
template void poisson::Physics::register_dirichlet<goal::Traits::Jacobian>(goal::FieldManager fm);
