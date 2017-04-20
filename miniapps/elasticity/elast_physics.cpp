#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_field.hpp>
#include <goal_indexer.hpp>
#include <goal_state_fields.hpp>
#include <Teuchos_ParameterList.hpp>

#include "elast_physics.hpp"

namespace elast {

using Teuchos::rcp;
using Teuchos::rcpFromRef;

static RCP<ParameterList> get_valid_params(RCP<goal::Discretization> d) {
  auto p = rcp(new ParameterList);
  p->sublist("dirichlet bcs");
  p->sublist("neumann bcs");
  p->sublist("qoi");
  for (int i = 0; i < d->get_num_elem_blocks(); ++i) {
    auto block = d->get_elem_block_name(i);
    p->sublist(block);
  }
  return p;
}

static void validate_params(
    RCP<const ParameterList> p, RCP<goal::Discretization> d) {
  assert(p->isSublist("dirichlet bcs"));
  for (int i = 0; i < d->get_num_elem_blocks(); ++i) {
    auto block = d->get_elem_block_name(i);
    assert(p->isSublist(block));
  }
  p->validateParameters(*get_valid_params(d), 0);
}

Physics::Physics(
    RCP<const ParameterList> p, RCP<goal::Discretization> d)
    : goal::Physics(d) {
  params = p;
  validate_params(p, disc);
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

RCP<const ParameterList> Physics::get_dbc_params() {
  return rcpFromRef(params->sublist("dirichlet bcs"));
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
  fm->writeGraphvizFile<Jacobian>("primal.dot", true, true);
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
  fm->writeGraphvizFile<Jacobian>("dual.dot", true, true);
}

void Physics::build_dual_neumann(FieldManager fm) {
  set_dual();
  register_neumann<Jacobian>(fm);
}

void Physics::build_dual_dirichlet(FieldManager fm) {
  set_dual();
  register_dirichlet<Jacobian>(fm);
}

void Physics::build_error_volumetric(FieldManager fm) {
  set_error();
  register_volumetric<Residual>(fm);
  fm->writeGraphvizFile<Residual>("error.dot", true, true);
}

void Physics::build_error_neumann(FieldManager fm) {
  set_error();
  register_neumann<Residual>(fm);
}

void Physics::build_error_dirichlet(FieldManager fm) {
  set_error();
  (void)fm;
}

}  /* namespace elast */

#include <goal_ev_utils.hpp>
#include <goal_ev_dirichlet_bcs.hpp>

#include "elast_ev_stress.hpp"
#include "elast_ev_von_mises.hpp"
#include "elast_ev_residual.hpp"
#include "elast_ev_traction_bcs.hpp"

template <typename EvalT>
void elast::Physics::register_volumetric(goal::FieldManager fm) {

  /* register the displacements. */
  RCP<goal::Field> disp = u[0];
  goal::register_dof<EvalT>(disp, indexer, fm);

  /* get the material model parameters. */
  auto eb_idx = indexer->get_elem_block_idx();
  auto eb_name = disc->get_elem_block_name(eb_idx);
  auto mat_params = rcpFromRef(params->sublist(eb_name));

  /* compute the Cauchy stress tensor. */
  auto cauchy =
    rcp(new elast::Stress<EvalT, goal::Traits>(disp, states, mat_params));
  fm->registerEvaluator<EvalT>(cauchy);

  /* compute the momentum residual. */
  if (is_primal || is_dual) {
    auto residual = rcp(new elast::Residual<EvalT, goal::Traits>(disp));
    fm->registerEvaluator<EvalT>(residual);
  }

  /* require the primal scatter. */
  if (is_primal)
    goal::require_primal_scatter<EvalT>(disp, indexer, fm);

  /* require the dual scatter. */
  if (is_dual) {
    goal::require_adjoint_scatter<EvalT>(disp, indexer, fm);
    auto vm = rcp(new elast::VonMises<EvalT, goal::Traits>(disp));
    fm->registerEvaluator<EvalT>(vm);
    auto qoip = rcpFromRef(params->sublist("qoi"));
    goal::require_qoi(qoip, disp, indexer, fm);
  }

  /* finalize the field manager registration. */
  goal::set_extended_data_type_dims(indexer, fm);
  fm->postRegistrationSetupForType<EvalT>(NULL);
}

template <typename EvalT>
void elast::Physics::register_neumann(goal::FieldManager fm) {
  /* bail if there are no neumann bcs specified. */
  if (!params->isSublist("neumann bcs")) return;

  /* build the traction bcs evaluator. */
  auto bc = rcpFromRef(params->sublist("neumann bcs"));
  auto ev = rcp(new elast::TractionBCs<EvalT, goal::Traits>(indexer, bc));
  fm->registerEvaluator<EvalT>(ev);
  fm->requireField<EvalT>(*ev->evaluatedFields()[0]);

  /* set the FAD data and finalize PHX field manager registration. */
  goal::set_extended_data_type_dims(indexer, fm);
  fm->postRegistrationSetupForType<EvalT>(NULL);
}

template <typename EvalT>
void elast::Physics::register_dirichlet(goal::FieldManager fm) {
  /* bail if we are estimating the error. */
  if (is_error) return;

  /* set up the Dirichlet BC evaluator. */
  auto dbc = get_dbc_params();
  auto ev = rcp(new goal::DirichletBCs<EvalT, goal::Traits>(
        dbc, indexer, true, is_dual));
  fm->registerEvaluator<EvalT>(ev);
  fm->requireField<EvalT>(*ev->evaluatedFields()[0]);

  /* set the FAD data and finalize PHX field manager registration. */
  goal::set_extended_data_type_dims(indexer, fm);
  fm->postRegistrationSetupForType<EvalT>(NULL);
}

template void elast::Physics::register_volumetric<goal::Traits::Residual>(goal::FieldManager fm);
template void elast::Physics::register_volumetric<goal::Traits::Jacobian>(goal::FieldManager fm);
template void elast::Physics::register_neumann<goal::Traits::Residual>(goal::FieldManager fm);
template void elast::Physics::register_neumann<goal::Traits::Jacobian>(goal::FieldManager fm);
template void elast::Physics::register_dirichlet<goal::Traits::Residual>(goal::FieldManager fm);
template void elast::Physics::register_dirichlet<goal::Traits::Jacobian>(goal::FieldManager fm);
