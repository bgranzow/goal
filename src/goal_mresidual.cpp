#include "goal_control.hpp"
#include "goal_displacement.hpp"
#include "goal_model.hpp"
#include "goal_mresidual.hpp"
#include "goal_vector_weight.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
MResidual<T>::MResidual(
  RCP<Integrator> u_,
  RCP<Integrator> w_, 
  RCP<Model<T>> model_) :
    u(rcp_static_cast<Displacement<T>>(u_)),
    w(rcp_static_cast<VectorWeight>(w_)),
    model(model_),
    num_dims(u->get_num_dims()) {
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(u));
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(w));
  this->name = "mresidual";
}

template <typename T>
void MResidual<T>::at_point(apf::Vector3 const&, double ipw, double dv) {
  auto P = model->get_first_pk();
  for (int n = 0; n < u->get_num_nodes(); ++n)
  for (int i = 0; i < num_dims; ++i)
  for (int j = 0; j < num_dims; ++j)
    u->resid(n, i) += P(i,j) * w->grad(n,i,j) * ipw * dv;
}

template class MResidual<ST>;
template class MResidual<FADT>;

}
