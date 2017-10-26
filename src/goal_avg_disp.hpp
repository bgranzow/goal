#ifndef goal_avg_disp_hpp
#define goal_avg_disp_hpp

#include "goal_qoi.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Integrator;
template <typename T> class Displacement;

template <typename T>
class AvgDisp : public QoI<T> {
  public:
    AvgDisp(RCP<Integrator> disp);
    void at_point(apf::Vector3 const&, double w, double dv);
  private:
    RCP<Displacement<T>> u;
    int num_dims;
};

}

#endif
