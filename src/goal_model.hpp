#ifndef goal_model_hpp
#define goal_model_hpp

#include <MiniTensor.h>
#include "goal_integrator.hpp"

namespace goal {

template <typename T>
class Model : public Integrator {
  private:
    using Tensor = minitensor::Tensor<T>;
  public:
    Model();
    virtual ~Model();
    virtual void set_elem_set(int) {}
    virtual void in_elem(apf::MeshElement*) {}
    virtual void at_point(apf::Vector3 const&, double, double) {}
    virtual void out_elem() {}
    virtual int get_num_dims() = 0;
    virtual Tensor& get_cauchy() = 0;
    virtual Tensor& get_first_pk() = 0;
};

}

#endif
