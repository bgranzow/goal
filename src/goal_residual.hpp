#ifndef goal_residual_hpp
#define goal_residual_hpp

#include "goal_integrator.hpp"

namespace goal {

class Disc;
class ScalarWeight;
template <typename T> class Soln;

template <typename T>
class Residual : public Integrator {
  public:
    Residual(
        RCP<Integrator> u,
        RCP<Integrator> w,
        std::string const& f);
    void pre_process(SolInfo* s);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const&, double, double);
    void out_elem();
    void post_process(SolInfo*);
  private:
    RCP<Soln<T>> u;
    RCP<ScalarWeight> w;
    std::string f;
    Disc* disc;
    apf::MeshElement* elem;
    int num_dims;
};

}

#endif
