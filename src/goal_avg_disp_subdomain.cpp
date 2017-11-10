#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_displacement.hpp"
#include "goal_avg_disp_subdomain.hpp"
#include "goal_sol_info.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("type", "");
  p.set<std::string>("elem set", "");
  return p;
}

template <typename T>
AvgDispSubdomain<T>::AvgDispSubdomain(
    ParameterList p, RCP<Integrator> disp) {
  params = p;
  params.validateParameters(get_valid_params(), 0);
  u = rcp_static_cast<Displacement<T>>(disp);
  num_dims = u->get_num_dims();
  es_idx = 0;
  this->name = "avg disp subdomain";
}

template <typename T>
void AvgDispSubdomain<T>::pre_process(SolInfo* s) {
  this->disc = s->get_disc();
  this->qoi_value = 0.0;
  auto esn = params.get<std::string>("elem set");
  es_idx = this->disc->get_elem_set_idx(esn);
}

template <typename T>
void AvgDispSubdomain<T>::set_elem_set(int idx) {
  if (idx == es_idx) op = &AvgDispSubdomain<T>::do_avg;
  else op = &AvgDispSubdomain<T>::do_null;
}

template <typename T>
void AvgDispSubdomain<T>::do_avg(double w, double dv) {
  for (int i = 0; i < num_dims; ++i)
    this->elem_value += u->val(i) * w * dv;
  this->elem_value /= num_dims;
}

template <typename T>
void AvgDispSubdomain<T>::do_null(double, double) {
}

template <typename T>
void AvgDispSubdomain<T>::at_point(
    apf::Vector3 const&, double w, double dv) {
  op(this, w, dv);
}

template class AvgDispSubdomain<ST>;
template class AvgDispSubdomain<FADT>;

}
