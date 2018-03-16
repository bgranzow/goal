#ifndef goal_elastic_hpp
#define goal_elastic_hpp

#include "goal_model.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class States;
template <typename T> class Displacement;

template <typename T>
class Elastic : public Model<T> {
  private:
    using Tensor = minitensor::Tensor<T>;
  public:
    Elastic(
        RCP<Integrator> disp,
        States* s,
        bool set,
        ParameterList const& p);
    void set_elem_set(int es_idx);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const&, double, double);
    void out_elem();
    int get_num_dims() { return num_dims; }
    Tensor& get_cauchy() { return sigma; }
    Tensor& get_first_pk() { return sigma; }
  private:
    RCP<Displacement<T>> u;
    States* states;
    ParameterList params;
    int num_dims;
    bool save_state;
    double mu;
    double lambda;
    Tensor eps, sigma, I;
    apf::MeshEntity* elem;
};

}

#endif
