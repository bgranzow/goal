#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_displacement.hpp"
#include "goal_kinematics.hpp"
#include "goal_model.hpp"
#include "goal_presidual.hpp"
#include "goal_pressure.hpp"
#include "goal_scalar_weight.hpp"
#include "goal_sol_info.hpp" 

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
PResidual<T>::PResidual(
  RCP<Integrator> disp,
  RCP<Integrator> pressure,
  RCP<Integrator> weight,
  RCP<Model<T>> model,
  RCP<Kinematics<T>> kinematics,
  ParameterList const& mat) :
    u(rcp_static_cast<Displacement<T>>(disp)),
    p(rcp_static_cast<Pressure<T>>(pressure)),
    w(rcp_static_cast<ScalarWeight>(weight)),
    m(model),
    k(kinematics),
    params(mat),
    disc(0),
    num_dims(p->get_num_dims()),
    kappa(0.0) {
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(p));
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(w));
  this->name = "pressure";
}

template <typename T>
void PResidual<T>::pre_process(SolInfo* s) {
  disc = s->get_disc();
}

template <typename T>
void PResidual<T>::set_elem_set(int es_idx) {
  auto esn = disc->get_elem_set_name(es_idx);
  ParameterList mat = params.sublist(esn);
  double E = mat.get<double>("E");
  double nu = mat.get<double>("nu");
  kappa = E / (3.0 * (1.0 - 2.0 * nu));
  mu = E / (2.0 * (1.0 + nu));
  lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

template <typename T>
void PResidual<T>::do_large_strain(double ipw, double dv) {
  auto J = k->get_det_def_grad();
  T dUdJ = 0.5*(J - 1.0/J);
  for (int n = 0; n < p->get_num_nodes(); ++n)
    p->resid(n) += ((p->val() / kappa) - dUdJ) * w->val(n) * ipw * dv;
}

template <typename T>
void PResidual<T>::do_small_strain(double ipw, double dv) {
  T pbar = 0.0;
  for (int i = 0; i < num_dims; ++i)
    pbar += u->grad(i, i);
  pbar *= (2.0*mu/num_dims + lambda);
  for (int n = 0; n < p->get_num_nodes(); ++n)
    p->resid(n) += (p->val() - pbar) * w->val(n) * ipw * dv;
}

template <typename T>
void PResidual<T>::at_point(apf::Vector3 const&, double ipw, double dv) {
  if (m->small_strain())
    do_small_strain(ipw, dv);
  else
    do_large_strain(ipw, dv);
}

template <typename T>
void PResidual<T>::post_process(SolInfo*) {
  disc = 0;
}

template class PResidual<ST>;
template class PResidual<FADT>;

}
