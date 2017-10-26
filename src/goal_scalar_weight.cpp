#include "goal_control.hpp"
#include "goal_scalar_weight.hpp"

namespace goal {

ScalarWeight::ScalarWeight(apf::Field* base) {
  GOAL_DEBUG_ASSERT(apf::getValueType(base) == apf::SCALAR);
  shape = apf::getShape(base);
  auto fname = (std::string)apf::getName(base);
  this->name = fname.substr(0, 1) + "w";
}

ScalarWeight::~ScalarWeight() {
}

ST const& ScalarWeight::val(const int node) const {
  return BF[node];
}

ST const& ScalarWeight::grad(const int node, const int i) const {
  return GBF[node][i];
}

void ScalarWeight::in_elem(apf::MeshElement* me) {
  elem = me;
}

void ScalarWeight::at_point(apf::Vector3 const& p, double, double) {
  apf::getBF(shape, elem, p, BF);
  apf::getGradBF(shape, elem, p, GBF);
}

void ScalarWeight::out_elem() {
  elem = 0;
}

}
