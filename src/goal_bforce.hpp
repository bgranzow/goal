#ifndef goal_bforce_hpp
#define goal_bforce_hpp

#include "goal_integrator.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class States;
class VectorWeight;

template <typename T>
class BForce : public Integrator {
  public:
    BForce(
        RCP<Integrator> u,
        RCP<Integrator> w,
        States* s,
        ParameterList const& p);
    void set_elem_set(int es_idx);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const& xi, double ipw, double dv);
    void out_elem();
  private:
    RCP<Displacement<T>> u;
    RCP<VectorWeight> w;
    States* states;
    ParameterList params;
    int num_dims;
    double mu;
    double k;
    apf::MeshElement* elem;
};

}

#endif
