#include <apf.h>
#include <apfMesh2.h>
#include <PCU.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_ks_vm.hpp"
#include "goal_model.hpp"
#include "goal_sol_info.hpp"
#include "goal_states.hpp"
#include "goal_von_mises.hpp"

namespace goal {

using Teuchos::rcp_static_cast;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("type", "");
  p.set<double>("rho", 0.0);
  return p;
}

template <typename T>
KSVM<T>::KSVM(ParameterList const& p_, RCP<Integrator> m_) {
  params = p_;
  params.validateParameters(get_valid_params(), 0);
  model = rcp_static_cast<Model<T>>(m_);
  num_dims = model->get_num_dims();
  rho = params.get<double>("rho");
  max = 0.0;
  scale = 0.0;
  this->name = "ks max vm";
}

static double get_max_vm(Disc* d) {
  double max_vm = 0.0;
  minitensor::Tensor<ST> sigma(3);
  sigma.fill(0.0);
  auto m = d->get_apf_mesh();
  auto dim = d->get_num_dims();
  apf::MeshEntity* elem;
  auto it = m->begin(dim);
  while ((elem = m->iterate(it))) {
    get_tensor(d, "sigma", elem, sigma);
    double vm = compute_von_mises(sigma);
    max_vm = std::max(max_vm, vm);
  }
  m->end(it);
  print(" > max vm: %.15e", max_vm);
  PCU_Max_Doubles(&max_vm, 1);
  GOAL_DEBUG_ASSERT(max_vm > 0);
  return max_vm;
}

static double get_scale(Disc* d, double max, double rho) {
  double scale = 0.0;
  minitensor::Tensor<ST> sigma(3);
  sigma.fill(0.0);
  apf::Vector3 xi;
  apf::MeshEntity* elem;
  auto m = d->get_apf_mesh();
  auto dim = d->get_num_dims();
  auto it = m->begin(dim);
  while ((elem = m->iterate(it))) {
    auto me = apf::createMeshElement(m, elem);
    apf::getIntPoint(me, 1, 0, xi);
    double dv = apf::getDV(me, xi);
    double w = apf::getIntWeight(me, 1, 0);
    get_tensor(d, "sigma", elem, sigma);
    double vm = compute_von_mises(sigma);
    scale += std::exp( rho * (vm - max) ) * w * dv;
    apf::destroyMeshElement(me);
  }
  m->end(it);
  PCU_Add_Doubles(&scale, 1);
  GOAL_DEBUG_ASSERT(scale > 0);
  return scale;
}

template <typename T>
void KSVM<T>::pre_process(SolInfo* s) {
  this->qoi_value = 0.0;
  this->disc = s->get_disc();
  max = get_max_vm(this->disc);
  scale = get_scale(this->disc, max, rho);
}

template <typename T>
void KSVM<T>::at_point(apf::Vector3 const&, double w, double dv) {
  minitensor::Tensor<T> sigma3x3(3);
  sigma3x3.fill(0.0);
  auto sigma = model->get_cauchy();
  for (int i = 0; i < num_dims; ++i)
  for (int j = 0; j < num_dims; ++j)
    sigma3x3(i,j) = sigma(i,j);
  T vm = compute_von_mises<T>(sigma3x3);
  this->elem_value += (1.0/(rho*scale)) * std::exp(rho * (vm - max)) * w * dv;
}

template <typename T>
void KSVM<T>::post_process(SolInfo*) {
  this->disc = 0;
  this->qoi_value = max + (1.0/rho)*std::log(scale);
}

template class KSVM<ST>;
template class KSVM<FADT>;

}
