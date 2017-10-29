#ifndef goal_J2_hpp
#define goal_J2_hpp

#include "goal_model.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class States;
template <typename T> class Kinematics;

template <typename T>
class J2 : public Model<T> {
  private:
    using Tensor = minitensor::Tensor<T>;
  public:
    J2(
        RCP<Kinematics<T>> k,
        States* s,
        bool save,
        ParameterList const& p);
    void set_elem_set(int es_idx);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const&, double, double);
    void out_elem();
    int get_num_dims() { return num_dims; }
    Tensor& get_cauchy() { return sigma; }
    Tensor& get_first_pk();
  private:
    RCP<Kinematics<T>> kinematics;
    States* states;
    ParameterList params;
    int num_dims;
    bool save_state;
    double E, nu, K, Y, kappa, mu, sq23;
    T Jm23, mubar, smag, f, dgam, eqps;
    Tensor Fp, Fpinv, Cpinv, be, s, N, Fpn, I;
    Tensor sigma, P;
    apf::MeshEntity* elem;
};

}

#endif
