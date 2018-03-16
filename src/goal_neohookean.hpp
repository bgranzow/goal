#ifndef goal_neohookean_hpp
#define goal_neohookean_hpp

#include "goal_model.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class States;
template <typename T> class Kinematics;

template <typename T>
class Neohookean : public Model<T> {
  private:
    using Tensor = minitensor::Tensor<T>;
  public:
    Neohookean(
        RCP<Kinematics<T>> kinematics,
        States* s,
        bool set,
        ParameterList const& p);
    void set_elem_set(int es_idx);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const&, double, double);
    void out_elem();
    bool small_strain() { return false; }
    int get_num_dims() { return num_dims; }
    Tensor& get_cauchy() { return sigma; }
    Tensor& get_first_pk();
  private:
    RCP<Kinematics<T>> k;
    States* states;
    ParameterList params;
    int num_dims;
    bool save_state;
    double mu;
    double kappa;
    Tensor b, I, sigma, P;
    apf::MeshEntity* elem;
};

}

#endif
