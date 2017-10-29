#include <apf.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_J2.hpp"
#include "goal_kinematics.hpp"
#include "goal_scalar_types.hpp"
#include "goal_states.hpp"

namespace goal {

static ParameterList get_valid_params(Disc* d) {
  ParameterList p;
  for (int es = 0; es < d->get_num_elem_sets(); ++es) {
    auto es_name = d->get_elem_set_name(es);
    p.sublist(es_name).set<double>("E", 0.0);
    p.sublist(es_name).set<double>("nu", 0.0);
    p.sublist(es_name).set<double>("K", 0.0);
    p.sublist(es_name).set<double>("Y", 0.0);
    p.sublist(es_name).set<double>("c0", 0.0);
  }
  return p;
}

template <typename T>
J2<T>::J2(
  RCP<Kinematics<T>> k,
  States* s,
  bool save,
  ParameterList const& p) :
    kinematics(k),
    states(s),
    params(p),
    num_dims(k->get_num_dims()),
    save_state(save),
    sq23(std::sqrt(2.0/3.0)),
    Fp(num_dims),
    Fpinv(num_dims),
    Cpinv(num_dims),
    be(num_dims),
    s(num_dims),
    N(num_dims),
    Fpn(num_dims),
    I(minitensor::eye<T>(num_dims)),
    sigma(num_dims),
    P(num_dims),
    elem(0) {
  auto d = states->get_disc();
  params.validateParameters(get_valid_params(d));
}

template <typename T>
void J2<T>::set_elem_set(int es_idx) {
  auto d = states->get_disc();
  auto es_name = d->get_elem_set_name(es_idx);
  ParameterList mat = params.sublist(es_name);
  E = mat.get<double>("E");
  nu = mat.get<double>("nu");
  K = mat.get<double>("K");
  Y = mat.get<double>("Y");
  kappa = E / (3.0 * (1.0 - 2.0 * nu));
  mu = E / (2.0 * (1.0 + nu));
}

template <typename T>
void J2<T>::in_elem(apf::MeshElement* me) {
  elem = apf::getMeshEntity(me);
}

template <typename T>
void J2<T>::at_point(apf::Vector3 const&, double, double) {

  // get the underlying discretization
  auto d = states->get_disc();

  // get def grad quantities
  auto F = kinematics->get_def_grad();
  auto J = kinematics->get_det_def_grad();
  Jm23 = std::pow(J, -2.0/3.0);

  // get plastic def grad quantities
  goal::get_tensor<T>(d, "Fp_old", elem, Fp);
  Fpinv = minitensor::inverse(Fp);

  // compute the trial state
  Cpinv = Fpinv * minitensor::transpose(Fpinv);
  be = Jm23 * F * Cpinv * minitensor::transpose(F);
  s = mu * minitensor::dev(be);
  mubar = minitensor::trace(be) * mu / num_dims;

  // check the yield condition
  T smag = minitensor::norm(s);
  goal::get_scalar<T>(d, "eqps_old", elem, eqps);
  f = smag - sq23 * (Y + K * eqps);

  // plastic increment
  if (f > 1.0e-12) {
    int iter = 0;
    bool converged = false;
    dgam = 0.0;
    T H(0.0), dH(0.0), alpha(0.0), res(0.0);

    T X = 0.0;
    T R = f;
    T dRdX = -2.0 * mubar * (1.0 + H / (3.0 * mubar));

    while ((! converged) && (iter < 30)) {
      iter++;
      X = X - R / dRdX;
      alpha = eqps + sq23 * X;
      H = K * alpha;
      dH = K;
      R = smag - (2.0 * mubar * X + sq23 * (Y + H));
      dRdX = -2.0 * mubar * (1.0 + dH / (3.0 * mubar));
      res = std::abs(R);
      if ((res < 1.0e-11) || (res / Y < 1.0e-11) || (res / f < 1.0e-11))
        converged = true;
      if (iter == 30)
        fail("J2: return mapping failed");
    }

    // updates
    dgam = X;
    N = (1.0 / smag) * s;
    s -= 2.0 * mubar * dgam * N;
    if (save_state) goal::set_scalar<T>(d, "eqps", elem, alpha);

    // get Fpn
    Fpn = minitensor::exp(dgam * N) * Fp;
    if (save_state) goal::set_tensor<T>(d, "Fp", elem, Fpn);
  }

  // otherwise elastic increment
  else
    if (save_state) goal::set_scalar<T>(d, "eqps", elem, eqps);

  // compute stress
  T p = 0.5*kappa*(J - 1.0/J);
  sigma = s/J + p*I;
  if (save_state) goal::set_tensor<T>(d, "sigma", elem, sigma);

}

template <typename T>
void J2<T>::out_elem() {
  elem = 0;
}

template <typename T>
minitensor::Tensor<T>& J2<T>::get_first_pk() {
  auto F = kinematics->get_def_grad();
  auto J = kinematics->get_det_def_grad();
  auto Finv = minitensor::inverse(F);
  P = J * sigma * minitensor::transpose(Finv);
  return P;
}

template class J2<ST>;
template class J2<FADT>;

}
