#ifndef goal_von_mises_hpp
#define goal_von_mises_hpp

#include <MiniTensor.h>

namespace goal {

template <typename T>
T compute_von_mises(minitensor::Tensor<T> const& sigma);

}

#endif
