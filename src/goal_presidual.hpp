#ifndef goal_presidual_hpp
#define goal_presidual_hpp

namespace goal {

class ScalarWeight;
template <typename T> class Kinematics;
template <typename T> class Pressure;

template <typename T>
class PResidual : public Integrator {
  public:
    PResidual(
        RCP<Integrator> p,
        RCP<Integrator> w,
        RCP<Kinematics<T>> k,
        ParameterList const& mat);
    void pre_process(SolInfo* s);
    void set_elem_set(int es_idx);
    void at_point(apf::Vector3 const&, double, double);
    void post_process(SolInfo*);
  private:
    RCP<Pressure<T>> p;
    RCP<ScalarWeight> w;
    RCP<Kinematics<T>> k;
    ParameterList params;
    Disc* disc;
    int num_dims;
    double kappa;
};

}

#endif
