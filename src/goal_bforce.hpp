#ifndef goal_bforce_hpp
#define goal_bforce_hpp

#include "goal_integrator.hpp"

namespace goal {

using Teuchos::RCP;

class VectorWeight;

template <typename T>
class BForce : public Integrator {
  public:
    BForce(RCP<Integrator> u, RCP<Integrator> w);
    void at_point(apf::Vector3 const& xi, double ipw, double dv);
  private:
    RCP<Displacement<T>> u;
    RCP<VectorWeight> w;
    int num_dims;
};

}

#endif
