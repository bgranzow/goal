#ifndef goal_vector_weight_hpp
#define goal_vector_weight_hpp

#include <apf.h>
#include "goal_integrator.hpp"
#include "goal_scalar_types.hpp"

namespace goal {

class VectorWeight : public Integrator {
  public:
    VectorWeight(apf::Field* base);
    virtual ST const& val(int node, int i) const;
    virtual ST const& grad(int node, int i, int j) const;
    virtual void in_elem(apf::MeshElement* me);
    virtual void at_point(apf::Vector3 const& p, double, double);
    virtual void out_elem();
  protected:
    apf::FieldShape* shape;
    apf::MeshElement* elem;
    apf::NewArray<ST> BF;
    apf::NewArray<apf::Vector3> GBF;
};

}

#endif
