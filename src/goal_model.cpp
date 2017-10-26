#include "goal_model.hpp"
#include "goal_scalar_types.hpp"

namespace goal {

template <typename T>
Model<T>::Model() {
  this->name = "model";
}

template <typename T>
Model<T>::~Model() {
}

template class Model<ST>;
template class Model<FADT>;

}
