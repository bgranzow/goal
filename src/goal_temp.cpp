#include <apf.h>
#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_kinematics.hpp"
#include "goal_model.hpp"
#include "goal_states.hpp"
#include "goal_temp.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("temp", "");
  p.set<double>("ref temp", 0.0);
  return p;
}

template <typename T>
Temp<T>::Temp(
  ParameterList const& tp,
  ParameterList const& mp,
  RCP<Model<T>> cm,
  RCP<Kinematics<T>> k,
  States* s,
  bool save) :
    kin(k),
    model(cm),
    states(s),
    save_state(save),
    alpha(0.0),
    ref_temp(0.0),
    time(0.0),
    I(minitensor::eye<T>(kin->get_num_dims())),
    elem(0),
    temp_params(tp),
    mat_params(mp) {
  temp_params.validateParameters(get_valid_params());
  ref_temp = temp_params.get<double>("ref temp");
  temp = temp_params.get<std::string>("temp");
  this->name = "temp";
}

template <typename T>
void Temp<T>::set_time(double t_now, double) {
  time = t_now;
}

template <typename T>
void Temp<T>::set_elem_set(int es_idx) {
  auto d = states->get_disc();
  auto es_name = d->get_elem_set_name(es_idx);
  ParameterList mat = mat_params.sublist(es_name);
  alpha = mat.get<double>("alpha");
}

template <typename T>
void Temp<T>::in_elem(apf::MeshElement* me) {
  elem = me;
}

template <typename ScalarT>
void Temp<ScalarT>::at_point(apf::Vector3 const& p, double, double) {
  apf::Vector3 x(0,0,0);
  apf::mapLocalToGlobal(elem, p, x);
  double T = eval(temp, x[0], x[1], x[2], time);
  auto& sigma = model->get_cauchy();
  auto J = kin->get_det_def_grad();
  sigma -= 3.0 * alpha * (1.0 + 1.0/(J*J)) * (T - ref_temp) * I;
  if (save_state) {
    auto d = states->get_disc();
    auto ent = apf::getMeshEntity(elem);
    goal::set_tensor<ScalarT>(d, "sigma", ent, sigma);
  }
}

template <typename T>
void Temp<T>::out_elem() {
  elem = 0;
}

template class Temp<ST>;
template class Temp<FADT>;

}
