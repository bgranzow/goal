#ifndef goal_bforce_hpp
#define goal_bforce_hpp

#include <functional>
#include "goal_integrator.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class States;
class VectorWeight;

template <typename T>
class BForce : public Integrator {
  public:
    BForce(
        RCP<Integrator> u,
        RCP<Integrator> w,
        std::string const& f,
        States* s,
        ParameterList const& p);
    void set_elem_set(int es_idx);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const& xi, double ipw, double dv);
    void out_elem();
  private:
    void eval_elastic_squared(apf::Vector3 const& x, apf::Vector3& b);
    std::function<void(BForce<T>*,
        apf::Vector3 const& X, apf::Vector3& b)> op;
    RCP<Displacement<T>> u;
    RCP<VectorWeight> w;
    std::string ftype;
    States* states;
    ParameterList params;
    int num_dims;
    double mu;
    double k;
    double lambda;
    apf::MeshElement* elem;
};

}

#endif
