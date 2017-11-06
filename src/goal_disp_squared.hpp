#ifndef goal_disp_squared_hpp
#define goal_disp_squared_hpp

#include "goal_qoi.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Integrator;
template <typename T> class Displacement;

template <typename T>
class DispSquared : public QoI<T> {
  public:
    DispSquared(RCP<Integrator> disp);
    void at_point(apf::Vector3 const&, double w, double dv);
  private:
    RCP<Displacement<T>> u;
    int num_dims;
};

}

#endif
