#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_ev_vector_shape.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename EVALT, typename TRAITS>
VectorShape<EVALT, TRAITS>::VectorShape(RCP<Field> f)
    : field(f),
      wdv(f->get_wdv_name(), f->get_scalar_ip_dl()),
      shape(f->get_basis_name(), f->get_weight_dl()),
      grad_shape(f->get_grad_basis_name(), f->get_grad_weight_dl()) {
  /* populate the index dimensions for this evaluator. */
  num_nodes = field->get_num_elem_nodes();
  num_ips = field->get_num_elem_ips();
  num_dims = field->get_num_dims();
  assert(field->get_value_type() == VECTOR);

  /* populate the dependency structure of this evaluator. */
  this->addEvaluatedField(wdv);
  this->addEvaluatedField(shape);
  this->addEvaluatedField(grad_shape);
  this->setName(field->get_basis_name());
}

template <typename EVALT, typename TRAITS>
void VectorShape<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  this->utils.setFieldData(wdv, fm);
  this->utils.setFieldData(shape, fm);
  this->utils.setFieldData(grad_shape, fm);
  (void)d;
}

template <typename EVALT, typename TRAITS>
void VectorShape<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  apf::Vector3 p;
  apf::NewArray<double> BF;
  apf::NewArray<apf::Vector3> GBF;
  auto q = field->get_q_degree();
  auto s = field->get_apf_basis();
  auto m = field->get_apf_mesh();
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto me = apf::createMeshElement(m, e);
    for (int ip = 0; ip < num_ips; ++ip) {
      apf::getIntPoint(me, field->get_q_degree(), ip, p);
      auto w = apf::getIntWeight(me, q, ip);
      wdv(elem, ip) = w * apf::getDV(me, p);
      apf::getBF(s, me, p, BF);
      apf::getGradBF(s, me, p, GBF);
      for (int node = 0; node < num_nodes; ++node) {
        for (int i = 0; i < num_dims; ++i) {
          shape(elem, node, ip, i) = BF[node];
          for (int j = 0; j < num_dims; ++j)
            grad_shape(elem, node, ip, i, j) = GBF[node][j];
        }
      }
    }
    apf::destroyMeshElement(me);
  }
}

template class VectorShape<goal::Traits::Residual, goal::Traits>;
template class VectorShape<goal::Traits::Jacobian, goal::Traits>;

}  // namespace goal
