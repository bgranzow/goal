#include <apf.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_displacement.hpp"
#include "goal_elastic.hpp"
#include "goal_scalar_types.hpp"
#include "goal_states.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

static ParameterList get_valid_params(Disc* d) {
  ParameterList p;
  for (int es = 0; es < d->get_num_elem_sets(); ++es) {
    auto es_name = d->get_elem_set_name(es);
    p.sublist(es_name).set<double>("E", 0.0);
    p.sublist(es_name).set<double>("nu", 0.0);
    p.sublist(es_name).set<double>("c0", 0.0);
    p.sublist(es_name).set<double>("alpha", 0.0);
  }
  return p;
}

template <typename T>
Elastic<T>::Elastic(
  RCP<Integrator> disp,
  States* s,
  bool save,
  ParameterList const& p) :
    u(rcp_static_cast<Displacement<T>>(disp)),
    states(s),
    params(p),
    num_dims(u->get_num_dims()),
    save_state(save),
    eps(num_dims),
    sigma(num_dims),
    I(minitensor::eye<T>(num_dims)),
    elem(0) {
  auto d = states->get_disc();
  params.validateParameters(get_valid_params(d));
}

template <typename T>
void Elastic<T>::set_elem_set(int es_idx) {
  auto d = states->get_disc();
  auto es_name = d->get_elem_set_name(es_idx);
  ParameterList mat = params.sublist(es_name);
  double E = mat.get<double>("E");
  double nu = mat.get<double>("nu");
  mu = E / (2.0 * (1.0 + nu));
  lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

template <typename T>
void Elastic<T>::in_elem(apf::MeshElement* me) {
  elem = apf::getMeshEntity(me);
}

template <typename T>
void Elastic<T>::at_point(apf::Vector3 const&, double, double) {
  auto d = states->get_disc();
  for (int i = 0; i < num_dims; ++i)
  for (int j = 0; j < num_dims; ++j)
    eps(i,j) = 0.5 * (u->grad(i,j) + u->grad(j,i));
  sigma = 2.0*mu*eps + lambda*minitensor::trace(eps)*I;
  if (save_state)
    goal::set_tensor<T>(d, "sigma", elem, sigma);
}

template <typename T>
void Elastic<T>::out_elem() {
  elem = 0;
}

template class Elastic<ST>;
template class Elastic<FADT>;

}
