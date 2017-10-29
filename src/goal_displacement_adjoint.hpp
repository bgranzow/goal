#ifndef goal_displacement_adjoint_hpp
#define goal_displacement_adjoint_hpp

#include "goal_vector_weight.hpp"

namespace goal {

class DisplacementAdjoint : public VectorWeight {
  public:
    DisplacementAdjoint(apf::Field* z);
    ST const& val(int node, int i) const;
    ST const& grad(int node, int i, int j) const;
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const& p, double, double);
    void out_elem();
  private:
    int num_dims;
    int num_nodes;
    apf::Field* field;
    apf::Element* z_elem;
    apf::MeshElement* elem;
    apf::NewArray<apf::Vector3> values;
    apf::NewArray<apf::Matrix3x3> gradients;
};

}

#endif
