#include "goal_control.hpp"
#include "goal_mixed.hpp"
#include "goal_model.hpp"
#include "goal_pressure.hpp"
#include "goal_states.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
Mixed<T>::Mixed(
  RCP<Integrator> pr,
  RCP<Model<T>> m,
  States* s,
  const bool save) :
    p(rcp_static_cast<Pressure<T>>(pr)),
    model(m),
    states(s),
    save_state(save),
    elem(0),
    num_dims(p->get_num_dims()) {
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(p));
  GOAL_DEBUG_ASSERT(Teuchos::nonnull(m));
  this->name = "mixed";
}

template <typename T>
void Mixed<T>::in_elem(apf::MeshElement* me) {
  elem = apf::getMeshEntity(me);
}

template <typename T>
void Mixed<T>::at_point(apf::Vector3 const&, double, double) {
  T pbar = 0.0;
  auto d = states->get_disc();
  auto& sigma = model->get_cauchy();
  auto pressure = p->val();
  for (int i = 0; i < num_dims; ++i)
    pbar += sigma(i, i);
  pbar /= num_dims;
  for (int i = 0; i < num_dims; ++i)
    sigma(i, i) += pressure - pbar;
  if (save_state)
    goal::set_tensor<T>(d, "sigma", elem, sigma);
}

template <typename T>
void Mixed<T>::out_elem() {
  elem = 0;
}

template class Mixed<ST>;
template class Mixed<FADT>;

}
