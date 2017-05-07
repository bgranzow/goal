#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_ev_dual_vector_weight.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;

template <typename EVALT, typename TRAITS>
DualVectorWeight<EVALT, TRAITS>::DualVectorWeight(
    RCP<Field> z_, RCP<Field> z_fine_) :
    z(z_),
    z_fine(z_fine_),
    w(z_fine->get_name(), z->get_PU_dl()),
    grad_w(z_fine->get_grad_name(), z->get_grad_PU_dl()) {
  /* make sure we're doing sane things. */
  assert(z->get_value_type() == VECTOR);
  assert(z_fine->get_value_type() == VECTOR);

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
void DualVectorWeight<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(w, fm);
  this->utils.setFieldData(grad_w, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void DualVectorWeight<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  apf::Vector3 p;
  apf::Vector3 z_val;
  apf::Vector3 z_fine_val;
  apf::Vector3 z_fine_min_z_val;
  apf::Matrix3x3 grad_z_val;
  apf::Matrix3x3 grad_z_fine_val;
  apf::Matrix3x3 grad_z_fine_min_z_val;
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
      apf::getVector(z_elem, p, z_val);
      apf::getVector(z_fine_elem, p, z_fine_val);
      apf::getVectorGrad(z_elem, p, grad_z_val);
      apf::getVectorGrad(z_fine_elem, p, grad_z_fine_val);
      z_fine_min_z_val = z_fine_val - z_val;
      grad_z_fine_min_z_val = grad_z_fine_val - grad_z_val;
      for (int vtx = 0; vtx < num_vtx; ++vtx) {
        for (int i = 0; i < num_dims; ++i) {
          w(elem, vtx, ip, i) = z_fine_min_z_val[i] * PU[vtx];
          for (int j = 0; j < num_dims; ++j)
            grad_w(elem, vtx, ip, i, j) =
              grad_z_fine_min_z_val[i][j] * PU[vtx] +
              z_fine_min_z_val[i] * GPU[vtx][j];
        }
      }
    }
    apf::destroyElement(z_fine_elem);
    apf::destroyElement(z_elem);
    apf::destroyMeshElement(me);
  }
}

template class DualVectorWeight<goal::Traits::Residual, goal::Traits>;

} /* namespace goal */
