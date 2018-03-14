#include "goal_control.hpp"
#include "goal_displacement.hpp"
#include "goal_vector_weight.hpp"
#include "goal_bforce.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
BForce<T>::BForce(
  RCP<Integrator> u_,
  RCP<Integrator> w_) :
    u(rcp_static_cast<Displacement<T>>(u_)),
    w(rcp_static_cast<VectorWeight>(w_)),
    num_dims(u->get_num_dims()) {
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(u));
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(w));
  this->name = "bforce";
}

template <typename T>
void BForce<T>::at_point(apf::Vector3 const& xi, double ipw, double dv) {
  (void)xi;
  (void)ipw;
  (void)dv;
}

template class BForce<ST>;
template class BForce<FADT>;

}
