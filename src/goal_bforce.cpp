#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_displacement.hpp"
#include "goal_states.hpp"
#include "goal_vector_weight.hpp"

#include "goal_bforce.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
BForce<T>::BForce(
  RCP<Integrator> u_,
  RCP<Integrator> w_,
  std::string const& f_,
  States* s_,
  ParameterList const& p_) :
    u(rcp_static_cast<Displacement<T>>(u_)),
    w(rcp_static_cast<VectorWeight>(w_)),
    ftype(f_),
    states(s_),
    params(p_),
    num_dims(u->get_num_dims()) {
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(u));
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(w));
  if (ftype == "elastic squared")
    op = &BForce<T>::eval_elastic_squared;
  else
    fail("unkown body force: %s", ftype.c_str());
  this->name = "bforce";
}

template <typename T>
void BForce<T>::set_elem_set(int es_idx) {
  auto d = states->get_disc();
  auto es_name = d->get_elem_set_name(es_idx);
  ParameterList mat = params.sublist(es_name);
  double E = mat.get<double>("E");
  double nu = mat.get<double>("nu");
  k = E / (3.0 * (1.0 - 2.0 * nu));
  mu = E / (2.0 * (1.0 + nu));
  lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

template <typename T>
void BForce<T>::in_elem(apf::MeshElement* me) {
  elem = me;
}

template <typename T>
void BForce<T>::eval_elastic_squared(
    apf::Vector3 const&, apf::Vector3& b) {
  b[0] = -4.0*mu - 2.0*lambda;
  b[1] = -4.0*mu - 2.0*lambda;
}

template <typename T>
void BForce<T>::at_point(apf::Vector3 const& xi, double ipw, double dv) {
  apf::Vector3 b(0,0,0);
  apf::Vector3 x(0,0,0);
  apf::mapLocalToGlobal(elem, xi, x);
  op(this, x, b);
  for (int n = 0; n < u->get_num_nodes(); ++n)
  for (int i = 0; i < num_dims; ++i)
    u->resid(n, i) -= b[i] * w->val(n, i) * ipw * dv;
}

template <typename T>
void BForce<T>::out_elem() {
  elem = 0;
}

template class BForce<ST>;
template class BForce<FADT>;

}
