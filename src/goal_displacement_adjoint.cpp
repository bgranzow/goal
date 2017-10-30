#include <apfMesh.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_displacement_adjoint.hpp"

namespace goal {

DisplacementAdjoint::DisplacementAdjoint(apf::Field* f) :
    VectorWeight(f) {
  field = f;
  GOAL_DEBUG_ASSERT(apf::getValueType(field) == apf::VECTOR);
  num_dims = apf::getMesh(field)->getDimension();
  num_nodes = 0;
  elem = 0;
}

ST const& DisplacementAdjoint::val(int node, int i) const {
  return values[node][i];
}

ST const& DisplacementAdjoint::grad(int node, int i, int j) const {
  return gradients[node][i][j];
}

void DisplacementAdjoint::in_elem(apf::MeshElement* me) {
  elem = me;
  z_elem = apf::createElement(field, me);
  auto m = apf::getMesh(field);
  auto ent = apf::getMeshEntity(elem);
  auto type = m->getType(ent);
  num_nodes = shape->getEntityShape(type)->countNodes();
  values.allocate(num_nodes);
  gradients.allocate(num_nodes);
}

void DisplacementAdjoint::at_point(apf::Vector3 const& p, double, double) {
  apf::Vector3 z;
  apf::Matrix3x3 grad_zT;
  apf::Matrix3x3 grad_z;
  apf::getBF(shape, elem, p, BF);
  apf::getGradBF(shape, elem, p, GBF);
  apf::getVector(z_elem, p, z);
  apf::getVectorGrad(z_elem, p, grad_zT);
  grad_z = apf::transpose(grad_zT);
  for (int n = 0; n < num_nodes; ++n) {
    for (int i = 0; i < num_dims; ++i) {
      values[n][i] = z[i]*BF[n];
      for (int j = 0; j < num_dims; ++j)
        gradients[n][i][j] = grad_z[i][j]*BF[n] + z[i]*GBF[n][j];
    }
  }
}

void DisplacementAdjoint::out_elem() {
  apf::destroyElement(z_elem);
  elem = 0;
}

}
