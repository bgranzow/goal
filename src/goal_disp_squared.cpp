#include "goal_control.hpp"
#include "goal_displacement.hpp"
#include "goal_disp_squared.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
DispSquared<T>::DispSquared(RCP<Integrator> disp) :
    u(rcp_static_cast<Displacement<T>>(disp)),
    num_dims(u->get_num_dims()) {
  this->name = "disp squared";
}

template <typename T>
void DispSquared<T>::at_point(apf::Vector3 const&, double w, double dv) {
  for (int i = 0; i < num_dims; ++i)
    this->elem_value += u->val(i) * u->val(i) * w * dv;
}

template class DispSquared<ST>;
template class DispSquared<FADT>;

}
