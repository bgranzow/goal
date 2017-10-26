#include <apf.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_kinematics.hpp"
#include "goal_neohookean.hpp"
#include "goal_scalar_types.hpp"
#include "goal_states.hpp"

namespace goal {

static ParameterList get_valid_params(Disc* d) {
  ParameterList p;
  for (int es = 0; es < d->get_num_elem_sets(); ++es) {
    auto es_name = d->get_elem_set_name(es);
    p.sublist(es_name).set<double>("E", 0.0);
    p.sublist(es_name).set<double>("nu", 0.0);
    p.sublist(es_name).set<double>("c0", 0.0);
  }
  return p;
}

template <typename T>
Neohookean<T>::Neohookean(
  RCP<Kinematics<T>> kinematics,
  States* s,
  const bool save,
  ParameterList const& p) :
    k(kinematics),
    states(s),
    params(p),
    num_dims(k->get_num_dims()),
    save_state(save),
    b(num_dims),
    I(minitensor::eye<T>(num_dims)),
    sigma(num_dims),
    elem(0) {
  auto d = states->get_disc();
  params.validateParameters(get_valid_params(d));
}

template <typename T>
void Neohookean<T>::set_elem_set(const int es_idx) {
  auto d = states->get_disc();
  auto es_name = d->get_elem_set_name(es_idx);
  ParameterList mat = params.sublist(es_name);
  double E = mat.get<double>("E");
  double nu = mat.get<double>("nu");
  kappa = E / (3.0 * (1.0 - 2.0 * nu));
  mu = E / (2.0 * (1.0 + nu));
}

template <typename T>
void Neohookean<T>::in_elem(apf::MeshElement* me) {
  elem = apf::getMeshEntity(me);
}

template <typename T>
void Neohookean<T>::at_point(apf::Vector3 const&, double, double) {
  auto d = states->get_disc();
  auto F = k->get_def_grad();
  auto J = k->get_det_def_grad();
  T Jm13 = 1.0 / std::cbrt(J);
  T Jm23 = Jm13 * Jm13;
  T Jm53 = Jm23 * Jm23 * Jm13;
  b = F*minitensor::transpose(F);
  T p = 0.5*kappa*(J - 1.0/J);
  sigma = mu*Jm53*minitensor::dev(b) + p*I;
  if (save_state)
    goal::set_tensor<T>(d, "sigma", elem, sigma);
}

template <typename T>
void Neohookean<T>::out_elem() {
  elem = 0;
}

template <typename T>
minitensor::Tensor<T>& Neohookean<T>::get_first_pk() {
  auto F = k->get_def_grad();
  auto J = k->get_det_def_grad();
  auto Finv = minitensor::inverse(F);
  P = J * sigma * minitensor::transpose(Finv);
  return P;
}

template class Neohookean<ST>;
template class Neohookean<FADT>;

}
