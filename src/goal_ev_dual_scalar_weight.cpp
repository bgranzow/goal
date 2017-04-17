#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_ev_dual_scalar_weight.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;

template <typename EVALT, typename TRAITS>
DualScalarWeight<EVALT, TRAITS>::DualScalarWeight(
    RCP<Field> z_, RCP<Field> z_fine_) :
    z(z_),
    z_fine(z_fine_),
    w(z->get_name(), z->get_PU_dl()),
    grad_w(z->get_grad_name(), z->get_grad_PU_dl()) {
  /* make sure we're doing sane stuff */
  assert(z->get_value_type() == SCALAR);
  assert(z_fine->get_value_type() == SCALAR);

  /* populate the index dimensions for this evaluator. */
  num_vtx = z_fine->get_num_elem_vtx();
  num_ips = z_fine->get_num_elem_ips();
  num_dims = z_fine->get_num_dims();

  /* populate the dependency structure of this evaluator. */
  this->addEvaluatedField(w);
  this->addEvaluatedField(grad_w);
  this->setName(z->get_name());
}

template <typename EVALT, typename TRAITS>
void DualScalarWeight<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(w, fm);
  this->utils.setFieldData(grad_w, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void DualScalarWeight<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  apf::Vector3 p;
  apf::Vector3 grad_z_val;
  apf::Vector3 grad_z_fine_val;
  apf::Vector3 grad_z_fine_min_z_val;
  apf::NewArray<double> PU;
  apf::NewArray<apf::Vector3> GPU;
  auto m = z->get_apf_mesh();
  auto apf_z = z->get_apf_field();
  auto apf_z_fine = z_fine->get_apf_field();
  auto s = apf::getLagrange(1);
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto me = apf::createMeshElement(m, e);
    auto z_elem = apf::createElement(apf_z, me);
    auto z_fine_elem = apf::createElement(apf_z_fine, me);
    for (int ip = 0; ip < num_ips; ++ip) {
      apf::getIntPoint(me, z_fine->get_q_degree(), ip, p);
      apf::getBF(s, me, p, PU);
      apf::getGradBF(s, me, p, GPU);
      auto z_val = apf::getScalar(z_elem, p);
      auto z_fine_val = apf::getScalar(z_fine_elem, p);
      auto z_fine_min_z_val = z_fine_val - z_val;
      apf::getGrad(z_elem, p, grad_z_val);
      apf::getGrad(z_fine_elem, p, grad_z_fine_val);
      grad_z_fine_min_z_val = grad_z_fine_val - grad_z_val;
      for (int vtx = 0; vtx < num_vtx; ++vtx) {
        w(elem, vtx, ip) = z_fine_min_z_val * PU[vtx];
        for (int dim = 0; dim < num_dims; ++dim)
          grad_w(elem, vtx, ip, dim) =
            grad_z_fine_min_z_val[dim] * PU[vtx] +
            z_fine_min_z_val * GPU[vtx][dim];
      }
    }
  }
}

template class DualScalarWeight<goal::Traits::Residual, goal::Traits>;

} /* namespace goal */

