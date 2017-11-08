#include <apf.h>
#include <PCU.h>

#include "goal_avg_vm.hpp"
#include "goal_control.hpp"
#include "goal_model.hpp"
#include "goal_von_mises.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("type", "");
  p.set<int>("elem set", 0);
  return p;
}

template <typename T>
AvgVM<T>::AvgVM(ParameterList const& p, RCP<Integrator> m) {
  params = p;
  params.validateParameters(get_valid_params(), 0);
  model = rcp_static_cast<Model<T>>(m);
  es_idx = params.get<int>("elem set");
  num_dims = model->get_num_dims();
  this->name = "avg vm";
}

template <typename T>
void AvgVM<T>::set_elem_set(int idx) {
  if (idx == es_idx) op = &AvgVM<T>::do_avg;
  else op = &AvgVM<T>::do_null;
}

template <typename T>
void AvgVM<T>::at_point(apf::Vector3 const&, double w, double dv) {
  op(this, w, dv);
}

template <typename T>
void AvgVM<T>::do_avg(double w, double dv) {
  minitensor::Tensor<T> sigma3x3(3);
  sigma3x3.fill(0.0);
  auto sigma = model->get_cauchy();
  for (int i = 0; i < num_dims; ++i)
  for (int j = 0; j < num_dims; ++j)
    sigma3x3(i, j) = sigma(i, j);
  T vm = compute_von_mises<T>(sigma3x3);
  this->elem_value += vm * w * dv;
}

template <typename T>
void AvgVM<T>::do_null(double, double) {
}

template class AvgVM<ST>;
template class AvgVM<FADT>;

}
