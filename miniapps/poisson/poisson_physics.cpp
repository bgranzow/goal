#include <goal_control.hpp>
#include <goal_field.hpp>
#include <goal_ev_utils.hpp>
#include <goal_ev_dirichlet_bcs.hpp>
#include <goal_ev_qoi_scalar_point.hpp>
#include <Teuchos_ParameterList.hpp>

#include "poisson_physics.hpp"
#include "poisson_ev_residual.hpp"

namespace poisson {

using Teuchos::rcp;

static RCP<ParameterList> get_valid_params() {
  auto p = rcp(new ParameterList);
  p->sublist("dirichlet bcs");
  p->set<std::string>("forcing function", "");
  p->set<std::string>("point set", "");
  return p;
}

static void validate_params(RCP<const ParameterList> p) {
  assert(p->isSublist("dirichlet bcs"));
  assert(p->isType<std::string>("forcing function"));
  assert(p->isType<std::string>("point set"));
  p->validateParameters(*get_valid_params(), 0);
}

Physics::Physics(
    RCP<const ParameterList> p, RCP<goal::Discretization> d)
    : goal::Physics(d) {
  params = p;
  validate_params(p);
  ff = params->get<std::string>("forcing function");
  set = params->get<std::string>("point set");
  goal::FieldInfo info =
    {disc, "u", 1, 1, goal::SCALAR, goal::HIERARCHIC};
  u.push_back(rcp(new goal::Field(&info)));
  u[0]->set_associated_dof_idx(0);
  u[0]->set_dof_status(true);
  u[0]->set_seed(1.0);
}

Physics::~Physics() {
}

RCP<const ParameterList> Physics::get_dbc_params() {
  return rcpFromRef(params->sublist("dirichlet bcs"));
}

void Physics::restrict_z_fine() {
  project_field(z[0], z_fine[0]);
}

void Physics::build_primal_volumetric(FieldManager fm) {
  auto uu = u[0];
  goal::register_dof<R>(uu, indexer, fm);
  goal::register_dof<J>(uu, indexer, fm);
  auto resid_R = rcp(new poisson::Residual<R, T>(uu, ff));
  auto resid_J = rcp(new poisson::Residual<J, T>(uu, ff));
  fm->registerEvaluator<R>(resid_R);
  fm->registerEvaluator<J>(resid_J);
  goal::require_primal_scatter<R>(uu, indexer, fm);
  goal::require_primal_scatter<J>(uu, indexer, fm);
  goal::set_extended_data_type_dims(indexer, fm);
  fm->postRegistrationSetupForType<R>(NULL);
  fm->postRegistrationSetupForType<J>(NULL);
  fm->writeGraphvizFile<R>("primal_volumetric", true, true);
}

void Physics::build_primal_dirichlet(FieldManager fm) {
  auto dbc = rcpFromRef(params->sublist("dirichlet bcs"));
  goal::require_dbc<R>(dbc, indexer, true, false, fm);
  goal::require_dbc<J>(dbc, indexer, true, false, fm);
  goal::set_extended_data_type_dims(indexer, fm);
  fm->postRegistrationSetupForType<R>(NULL);
  fm->postRegistrationSetupForType<J>(NULL);
}

void Physics::build_dual_volumetric(FieldManager fm) {
  auto uu = u_fine[0];
  goal::register_dof<J>(uu, indexer, fm);
  auto resid = rcp(new poisson::Residual<J, T>(uu, ff));
  fm->registerEvaluator<J>(resid);
  goal::require_adjoint_scatter<J>(uu, indexer, fm);
  goal::require_qoi_scalar_point(uu, indexer, set, fm);
  goal::set_extended_data_type_dims(indexer, fm);
  fm->postRegistrationSetupForType<J>(NULL);
  fm->writeGraphvizFile<J>("dual_volumetric", true, true);
}

void Physics::build_dual_dirichlet(FieldManager fm) {
  auto dbc = rcpFromRef(params->sublist("dirichlet bcs"));
  goal::require_dbc<J>(dbc, indexer, true, true, fm);
  goal::set_extended_data_type_dims(indexer, fm);
  fm->postRegistrationSetupForType<J>(NULL);
}

void Physics::build_error_volumetric(FieldManager fm) {
  auto ee = e[0];
  auto uu = u_fine[0];
  auto zz = z_fine[0];
  auto zz_coarse = z[0];
  goal::register_dof<R>(uu, indexer, fm);
  goal::register_dual<R>(zz_coarse, zz, fm);
  auto resid = rcp(new poisson::Residual<R, T>(uu, zz, ff));
  fm->registerEvaluator<R>(resid);
  goal::require_error(uu, ee, indexer, fm);
  goal::set_extended_data_type_dims(indexer, fm);
  fm->postRegistrationSetupForType<R>(NULL);
  fm->writeGraphvizFile<R>("error_volumetric", true, true);
}

}  /* namespace poisson */
