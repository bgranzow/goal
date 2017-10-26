#include "goal_displacement.hpp"
#include "goal_kinematics.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

template <typename T>
Kinematics<T>::Kinematics(RCP<Integrator> disp) :
    u(rcp_static_cast<Displacement<T>>(disp)),
    num_dims(u->get_num_dims()),
    J(0.0),
    F(num_dims) {
  this->name = "kinematics";
}

template <typename T>
void Kinematics<T>::at_point(apf::Vector3 const&, double, double) {
  for (int i = 0; i < num_dims; ++i)
  for (int j = 0; j < num_dims; ++j)
    F(i, j) = u->grad(i, j);
  for (int i = 0; i < num_dims; ++i)
    F(i, i) += 1.0;
  J = minitensor::det(F);
}

template class Kinematics<ST>;
template class Kinematics<FADT>;

}
