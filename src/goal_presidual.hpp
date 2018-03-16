#ifndef goal_presidual_hpp
#define goal_presidual_hpp

namespace goal {

class ScalarWeight;
template <typename T> class Displacement;
template <typename T> class Kinematics;
template <typename T> class Model;
template <typename T> class Pressure;

template <typename T>
class PResidual : public Integrator {
  public:
    PResidual(
        RCP<Integrator> u,
        RCP<Integrator> p,
        RCP<Integrator> w,
        RCP<Model<T>> cm,
        RCP<Kinematics<T>> k,
        ParameterList const& mat);
    void pre_process(SolInfo* s);
    void set_elem_set(int es_idx);
    void at_point(apf::Vector3 const&, double, double);
    void post_process(SolInfo*);
  private:
    void do_small_strain(double ipw, double dv);
    void do_large_strain(double ipw, double dv);
    RCP<Displacement<T>> u;
    RCP<Pressure<T>> p;
    RCP<ScalarWeight> w;
    RCP<Model<T>> m;
    RCP<Kinematics<T>> k;
    ParameterList params;
    Disc* disc;
    int num_dims;
    double kappa, mu, lambda;
};

}

#endif
