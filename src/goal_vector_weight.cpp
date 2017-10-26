#include "goal_control.hpp"
#include "goal_vector_weight.hpp"

namespace goal {

VectorWeight::VectorWeight(apf::Field* base) {
  GOAL_DEBUG_ASSERT(apf::getValueType(base) == apf::VECTOR);
  shape = apf::getShape(base);
  auto fname = (std::string)apf::getName(base);
  this->name = fname.substr(0, 1) + "w";
}

ST const& VectorWeight::val(const int node, const int) const {
  return BF[node];
}

ST const& VectorWeight::grad(const int node, const int, const int j) const {
  return GBF[node][j];
}

void VectorWeight::in_elem(apf::MeshElement* me) {
  elem = me;
}

void VectorWeight::at_point(apf::Vector3 const& p, double, double) {
  apf::getBF(shape, elem, p, BF);
  apf::getGradBF(shape, elem, p, GBF);
}

void VectorWeight::out_elem() {
  elem = 0;
}

}
