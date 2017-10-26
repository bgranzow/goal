#ifndef goal_mresidual_hpp
#define goal_mresidual_hpp

#include "goal_integrator.hpp"

namespace goal {

using Teuchos::RCP;

class VectorWeight;
template <typename T> class Model;
template <typename T> class Displacement;

template <typename T>
class MResidual : public Integrator {
  public:
    MResidual(RCP<Integrator> u, RCP<Integrator> w, RCP<Model<T>> m);
    void at_point(apf::Vector3 const&, double ipw, double dv);
  private:
    RCP<Displacement<T>> u;
    RCP<VectorWeight> w;
    RCP<Model<T>> model;
    int num_dims;
};

}

#endif
