#include "goal_avg_disp.hpp"
#include "goal_control.hpp"
#include "goal_displacement.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
AvgDisp<T>::AvgDisp(RCP<Integrator> disp) :
    u(rcp_static_cast<Displacement<T>>(disp)),
    num_dims(u->get_num_dims()) {
  this->name = "avg disp";
}

template <typename T>
void AvgDisp<T>::at_point(apf::Vector3 const&, double w, double dv) {
  for (int i = 0; i < num_dims; ++i)
    this->elem_value += u->val(i) * w * dv;
  this->elem_value /= num_dims;
}

template class AvgDisp<ST>;
template class AvgDisp<FADT>;

}
