#include <apf.h>
#include <apfMesh.h>
#include <goal_control.hpp>
#include <goal_field.hpp>
#include <goal_traits.hpp>
#include <goal_workset.hpp>

#include "poisson_ev_residual.hpp"

namespace poisson {

template <typename EVALT, typename TRAITS>
Residual<EVALT, TRAITS>::Residual(
    RCP<goal::Field> f, std::string const& source) :
      sol(f),
      ff(source),
      wdv(f->get_wdv_name(), f->get_scalar_ip_dl()),
      w(f->get_basis_name(), f->get_weight_dl()),
      grad_w(f->get_grad_basis_name(), f->get_grad_weight_dl()),
      grad_u(f->get_grad_name(), f->get_grad_interpolated_dl()),
      resid(f->get_residual_name(), f->get_dl()) {
  /* populate the index dimensions for this evaluator. */
  num_nodes = sol->get_num_elem_nodes();
  num_ips = sol->get_num_elem_ips();
  num_dims = sol->get_num_dims();

  /* populate the dependency structure of this evaluator. */
  this->addDependentField(wdv);
  this->addDependentField(w);
  this->addDependentField(grad_w);
  this->addDependentField(grad_u);
  this->addEvaluatedField(resid);
  this->setName("Poisson Residual");
}

template <typename EVALT, typename TRAITS>
Residual<EVALT, TRAITS>::Residual(
    RCP<goal::Field> u, RCP<goal::Field> z, std::string const& source)
    : sol(u),
      ff(source),
      wdv(u->get_wdv_name(), u->get_scalar_ip_dl()),
      w(z->get_name(), z->get_PU_dl()),
      grad_w(z->get_grad_name(), z->get_grad_PU_dl()),
      grad_u(u->get_grad_name(), u->get_grad_interpolated_dl()),
      resid(u->get_residual_name(), z->get_residual_PU_dl()) {
  /* populate the index dimensions. */
  num_nodes = sol->get_num_elem_vtx();
  num_ips = sol->get_num_elem_ips();
  num_dims = sol->get_num_dims();

  /* populate the dependency structure of this evaluator. */
  this->addDependentField(wdv);
  this->addDependentField(w);
  this->addDependentField(grad_w);
  this->addDependentField(grad_u);
  this->addEvaluatedField(resid);
  this->setName("Poisson Residual");
}

template <typename EVALT, typename TRAITS>
void Residual<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(wdv, fm);
  this->utils.setFieldData(w, fm);
  this->utils.setFieldData(grad_w, fm);
  this->utils.setFieldData(grad_u, fm);
  this->utils.setFieldData(resid, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void Residual<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  apf::Vector3 p(0, 0, 0);
  apf::Vector3 x(0, 0, 0);
  auto m = sol->get_apf_mesh();
  auto q = sol->get_q_degree();
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto me = apf::createMeshElement(m, e);
    for (int node = 0; node < num_nodes; ++node)
      resid(elem, node) = ScalarT(0.0);
    for (int ip = 0; ip < num_ips; ++ip)
      for (int node = 0; node < num_nodes; ++node)
        for (int dim = 0; dim < num_dims; ++dim)
          resid(elem, node) += grad_u(elem, ip, dim) *
                               grad_w(elem, node, ip, dim) * wdv(elem, ip);
    for (int node = 0; node < num_nodes; ++node) {
      for (int ip = 0; ip < num_ips; ++ip) {
        apf::getIntPoint(me, q, ip, p);
        apf::mapLocalToGlobal(me, p, x);
        auto f = goal::eval(ff, x[0], x[1], x[2], 0.0);
        resid(elem, node) -= f * w(elem, node, ip) * wdv(elem, ip);
      }
    }
    apf::destroyMeshElement(me);
  }
}

template class Residual<goal::Traits::Residual, goal::Traits>;
template class Residual<goal::Traits::Jacobian, goal::Traits>;

} /* namespace poisson */
