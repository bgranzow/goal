#include "goal_von_mises.hpp"
#include "goal_scalar_types.hpp"

namespace goal {

template <typename T>
T compute_von_mises(minitensor::Tensor<T> const& sigma)
{
  T s1 = (sigma(0,0)-sigma(1,1))*(sigma(0,0)-sigma(1,1));
  T s2 = (sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2));
  T s3 = (sigma(2,2)-sigma(0,0))*(sigma(2,2)-sigma(0,0));
  T s4 = sigma(0,1)*sigma(0,1);
  T s5 = sigma(1,2)*sigma(1,2);
  T s6 = sigma(2,0)*sigma(2,0);
  T s7 = 0.5*(s1+s2+s3+6.0*(s4+s5+s6));
  return std::sqrt(s7);
}

template ST compute_von_mises(minitensor::Tensor<ST> const& s);
template FADT compute_von_mises(minitensor::Tensor<FADT> const& s);

}
