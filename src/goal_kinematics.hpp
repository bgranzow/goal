#ifndef goal_kinematics_hpp
#define goal_kinematics_hpp

#include <MiniTensor.h>
#include <Teuchos_RCP.hpp>
#include "goal_integrator.hpp"

namespace goal {

using Teuchos::RCP;

template <typename T> class Displacement;

template <typename T>
class Kinematics : public Integrator {
  private:
    using Tensor = minitensor::Tensor<T>;
  public:
    Kinematics(RCP<Integrator> disp);
    int get_num_dims() const { return num_dims; }
    T const& get_det_def_grad() const { return J; }
    Tensor const& get_def_grad() const { return F; }
    void at_point(apf::Vector3 const&, double, double);
  private:
    RCP<Displacement<T>> u;
    int num_dims;
    T J;
    Tensor F;
};

}

#endif
