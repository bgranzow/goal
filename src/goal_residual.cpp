#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_residual.hpp"
#include "goal_soln.hpp"
#include "goal_sol_info.hpp" 
#include "goal_scalar_weight.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
Residual<T>::Residual(
      RCP<Integrator> u_,
      RCP<Integrator> w_,
      std::string const& f_) :
    u(rcp_static_cast<Soln<T>>(u_)),
    w(rcp_static_cast<ScalarWeight>(w_)),
    f(f_),
    disc(0),
    elem(0),
    num_dims(u->get_num_dims()) {
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(u));
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(w));
  this->name = "residual";
}

template <typename T>
void Residual<T>::pre_process(SolInfo* s) {
  disc = s->get_disc();
}

template <typename T>
void Residual<T>::in_elem(apf::MeshElement* me) {
  elem = me;
}

template <typename T>
void Residual<T>::at_point(apf::Vector3 const& p, double ipw, double dv) {
  apf::Vector3 x(0,0,0);
  apf::mapLocalToGlobal(elem, p, x);
  double fval = eval(f, x[0], x[1], x[2], 0.0);
  for (int n = 0; n < u->get_num_nodes(); ++n)
    for (int i = 0; i < num_dims; ++i)
      u->resid(n) += u->grad(i) * w->grad(n, i) * ipw * dv;
  for (int n = 0; n < u->get_num_nodes(); ++n)
    u->resid(n) -= fval * w->val(n) * ipw * dv;


  for (int n = 0; n < u->get_num_nodes(); ++n)
    std::cout << u->resid(n) << std::endl;
}

template <typename T>
void Residual<T>::out_elem() {
  elem = 0;
}

template <typename T>
void Residual<T>::post_process(SolInfo*) {
  disc = 0;
}

template class Residual<ST>;
template class Residual<FADT>;

}
