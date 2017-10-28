#ifndef goal_pressure_adjoint_hpp
#define goal_pressure_adjoint_hpp

#include "goal_scalar_weight.hpp"

namespace goal {

class PressureAdjoint : public ScalarWeight {
  public:
    PressureAdjoint(apf::Field* z);
    ST const& val(int node) const;
    ST const& grad(int node, int i) const;
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const& p, double, double);
    void out_elem();
  private:
    int num_dims;
    int num_nodes;
    apf::Field* field;
    apf::Element* z_elem;
    apf::MeshElement* elem;
    apf::NewArray<ST> values;
    apf::NewArray<apf::Vector3> gradients;
};

}

#endif
