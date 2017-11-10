#ifndef goal_avg_disp_subdomain_hpp
#define goal_avg_disp_subdomain_hpp

#include "goal_qoi.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Integrator;
template <typename T> class Displacement;

template <typename T>
class AvgDispSubdomain : public QoI<T> {
  public:
    AvgDispSubdomain(ParameterList p, RCP<Integrator> disp);
    void pre_process(SolInfo* s);
    void set_elem_set(int es_idx);
    void at_point(apf::Vector3 const&, double w, double dv);
  private:
    void do_avg(double w, double dv);
    void do_null(double w, double dv);
    std::function<void(AvgDispSubdomain<T>*, double, double)> op;
    ParameterList params;
    RCP<Displacement<T>> u;
    int num_dims;
    int es_idx;
};

}

#endif
