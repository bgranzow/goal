#include <apfMesh.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_coarse_pressure_adjoint.hpp"

namespace goal {

CoarsePressureAdjoint::CoarsePressureAdjoint(apf::Field* f) :
   ScalarWeight(f) {
  field = f;
  GOAL_DEBUG_ASSERT(apf::getValueType(field) == apf::SCALAR);
  num_dims = 0;
  elem = 0;
  auto fname = (std::string)apf::getName(field);
  this->name = fname.substr(0,1) + "wc";
}

ST const& CoarsePressureAdjoint::val(int node) const {
  return values[node];
}

ST const& CoarsePressureAdjoint::grad(int node, int i) const {
  return gradients[node][i];
}

void CoarsePressureAdjoint::in_elem(apf::MeshElement* me) {
  elem = me;
  z_elem = apf::createElement(field, me);
  auto m = apf::getMesh(field);
  auto ent = apf::getMeshEntity(elem);
  auto type = m->getType(ent);
  num_nodes = shape->getEntityShape(type)->countNodes();
  values.allocate(num_nodes);
  gradients.allocate(num_nodes);
}

void CoarsePressureAdjoint::at_point(apf::Vector3 const& p, double, double) {
  apf::Vector3 grad_z;
  apf::getBF(shape, elem, p, BF);
  apf::getGradBF(shape, elem, p, GBF);
  double z = apf::getScalar(z_elem, p);
  apf::getGrad(z_elem, p, grad_z);
  for (int n = 0; n < num_nodes; ++n) {
    values[n] = z * BF[n];
    for (int i = 0; i < num_dims; ++i)
      gradients[n][i] = grad_z[i]*BF[n] + z*GBF[n][i];
  }
}

void CoarsePressureAdjoint::out_elem() {
  elem = 0;
  apf::destroyElement(z_elem);
}

}
