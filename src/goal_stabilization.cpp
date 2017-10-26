#include <apf.h>
#include <apfMesh2.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_kinematics.hpp"
#include "goal_model.hpp"
#include "goal_pressure.hpp"
#include "goal_scalar_weight.hpp"
#include "goal_sol_info.hpp"
#include "goal_stabilization.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
Stabilization<T>::Stabilization(
  RCP<Integrator> pressure,
  RCP<Integrator> weight,
  RCP<Model<T>> model,
  RCP<Kinematics<T>> kinematics,
  ParameterList const& mat) :
    p(rcp_static_cast<Pressure<T>>(pressure)),
    w(rcp_static_cast<ScalarWeight>(weight)),
    m(model),
    k(kinematics),
    params(mat),
    disc(0),
    mesh(0),
    elem(0),
    num_dims(p->get_num_dims()),
    mu(0.0),
    c0(0.0) {
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(p));
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(w));
  this->name = "stabilization";
}

template <typename T>
void Stabilization<T>::pre_process(SolInfo* s) {
  disc = s->get_disc();
  mesh = disc->get_apf_mesh();
}

template <typename T>
void Stabilization<T>::set_elem_set(const int es_idx) {
  auto es_name = disc->get_elem_set_name(es_idx);
  ParameterList mat = params.sublist(es_name);
  double E = mat.get<double>("E");
  double nu = mat.get<double>("nu");
  c0 = mat.get<double>("c0");
  mu = E / (2.0 * (1.0 + nu));
}

template <typename T>
void Stabilization<T>::in_elem(apf::MeshElement* me) {
  elem = apf::getMeshEntity(me);
}

static double get_size(apf::Mesh* m, apf::MeshEntity* e) {
  double h = 0.0;
  apf::Downward edges;
  int ne = m->getDownward(e, 1, edges);
  for (int i = 0; i < ne; ++i)
    h += apf::measure(m, edges[i]) * apf::measure(m, edges[i]);
  return std::sqrt(h/ne);
}

template <typename T>
void Stabilization<T>::at_point(apf::Vector3 const&, double ipw, double dv) {
  using minitensor::inverse;
  using minitensor::transpose;
  double h = get_size(mesh, elem);
  double tau = 0.5*c0*h*h/mu;
  auto J = k->get_det_def_grad();
  auto F = k->get_def_grad();
  auto Cinv = inverse(transpose(F)*F);
  for (int n = 0; n < p->get_num_nodes(); ++n)
  for (int i = 0; i < num_dims; ++i)
  for (int j = 0; j < num_dims; ++j)
    p->resid(n) -= tau * J * Cinv(i, j) * p->grad(i) * w->grad(n, j) * ipw * dv;
}

template <typename T>
void Stabilization<T>::out_elem() {
  elem = 0;
}

template <typename T>
void Stabilization<T>::post_process(SolInfo*) {
  disc = 0;
  mesh = 0;
}

template class Stabilization<ST>;
template class Stabilization<FADT>;

}
