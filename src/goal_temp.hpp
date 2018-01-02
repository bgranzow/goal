#ifndef goal_temp_hpp
#define goal_temp_hpp

#include <Teuchos_RCP.hpp>
#include "goal_integrator.hpp"

namespace goal {

using Teuchos::ParameterList;

class States;
template <typename T> class Model;
template <typename T> class Kinematics;

template <typename T>
class Temp : public Integrator {
  public:
    Temp(
        ParameterList const& tp,
        ParameterList const& mp,
        RCP<Model<T>> cm,
        RCP<Kinematics<T>> k,
        States* s,
        bool save);
    void set_time(double t_now, double t_old);
    void set_elem_set(int es_idx);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const&, double, double);
    void out_elem();
  private:
    RCP<Kinematics<T>> kin;
    RCP<Model<T>> model;
    States* states;
    bool save_state;
    double alpha;
    double ref_temp;
    double time;
    std::string temp;
    minitensor::Tensor<T> I;
    apf::MeshElement* elem;
    ParameterList temp_params;
    ParameterList mat_params;
};

}

#endif
