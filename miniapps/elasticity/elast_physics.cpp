#include <goal_control.hpp>
#include <goal_field.hpp>
#include <goal_state_fields.hpp>
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
  assert(p->isSublist("dirichlet bcs"));
  p->validateParameters(*get_valid_params(), 0);
}

Physics::Physics(
    RCP<const ParameterList> p, RCP<goal::Discretization> d)
    : goal::Physics(d) {
  params = p;
  validate_params(p);
  is_primal = false;
  is_dual = false;
  is_error = false;
  q_degree = 2;
  p_order = 2;
  goal::FieldInfo i =
    {disc, "u", q_degree, p_order, goal::VECTOR, goal::LAGRANGE};
  u.push_back(rcp(new goal::Field(&i)));
  u.back()->set_associated_dof_idx(0);
  u.back()->set_dof_status(true);
  u.back()->set_seed(1.0);
  states = rcp(new goal::StateFields(disc, q_degree));
  states->add("cauchy", goal::TENSOR, false);
}

Physics::~Physics() {
  for (std::size_t i = 0; i < u.size(); ++i) u[i] = Teuchos::null;
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

void Physics::build_primal_neumann(FieldManager fm) {
  set_primal();
  register_neumann<Residual>(fm);
  register_neumann<Jacobian>(fm);
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

void Physics::build_dual_neumann(FieldManager fm) {
  set_dual();
  (void)fm;
}

void Physics::build_dual_dirichlet(FieldManager fm) {
  set_dual();
  register_dirichlet<Jacobian>(fm);
}

void Physics::build_error_volumetric(FieldManager fm) {
  set_error();
  register_volumetric<Residual>(fm);
}

void Physics::build_error_neumann(FieldManager fm) {
  set_error();
  register_neumann<Residual>(fm);
}

void Physics::build_error_dirichlet(FieldManager fm) {
  set_error();
  (void)fm;
}

}  // namespace elast

#include <goal_ev_utils.hpp>

template <typename EvalT>
void elast::Physics::register_volumetric(goal::FieldManager fm) {
  (void)fm;
}

template <typename EvalT>
void elast::Physics::register_neumann(goal::FieldManager fm) {
  (void)fm;
}

template <typename EvalT>
void elast::Physics::register_dirichlet(goal::FieldManager fm) {
  (void)fm;
}

template void elast::Physics::register_volumetric<goal::Traits::Residual>(goal::FieldManager fm);
template void elast::Physics::register_volumetric<goal::Traits::Jacobian>(goal::FieldManager fm);
template void elast::Physics::register_neumann<goal::Traits::Residual>(goal::FieldManager fm);
template void elast::Physics::register_neumann<goal::Traits::Jacobian>(goal::FieldManager fm);
template void elast::Physics::register_dirichlet<goal::Traits::Residual>(goal::FieldManager fm);
template void elast::Physics::register_dirichlet<goal::Traits::Jacobian>(goal::FieldManager fm);
