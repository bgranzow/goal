#ifndef goal_stabilization_hpp
#define goal_stabilization_hpp

#include "goal_integrator.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Disc;
class Integrator;
class ScalarWeight;
template <typename T> class Kinematics;
template <typename T> class Model;
template <typename T> class Pressure;

template <typename T>
class Stabilization : public Integrator {
  public:
    Stabilization(
        RCP<Integrator> p,
        RCP<Integrator> w,
        RCP<Model<T>> m,
        RCP<Kinematics<T>> k,
        ParameterList const& mat);
    void pre_process(SolInfo*);
    void set_elem_set(int es_idx);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const&, double ipw, double dv);
    void out_elem();
    void post_process(SolInfo*);
  private:
    RCP<Pressure<T>> p;
    RCP<ScalarWeight> w;
    RCP<Model<T>> m;
    RCP<Kinematics<T>> k;
    ParameterList params;
    Disc* disc;
    apf::Mesh* mesh;
    apf::MeshEntity* elem;
    int num_dims;
    double mu;
    double c0;
};

}

#endif
