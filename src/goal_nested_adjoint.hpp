#ifndef goal_nested_adjoint_hpp
#define goal_nested_adjoint_hpp

#include <Teuchos_ParameterList.hpp>

namespace goal {

class Disc;
class Integrator;
class Mechanics;
class Nested;
class Primal;
class SolInfo;

using Teuchos::RCP;
using Teuchos::ParameterList;
using Evaluators = std::vector<RCP<Integrator>>;

class NestedAdjoint {
  public:
    NestedAdjoint(ParameterList const& p, Primal* pr);
    ~NestedAdjoint();
    void build_data();
    void run(double t_now, double t_old);
    void destroy_data();
  private:
    void print_banner(double t_now);
    void compute_adjoint(double t_now, double t_old);
    void solve(double t_now, double t_old);
    ParameterList params;
    Primal* primal;
    Disc* base_disc;
    Nested* nested_disc;
    Mechanics* mech;
    SolInfo* sol_info;
    Evaluators adjoint;
    apf::Field* z_displacement;
    apf::Field* z_pressure;
};

NestedAdjoint* create_nested_adjoint(ParameterList const& p, Primal* pr);
void destroy_nested_adjoint(NestedAdjoint* a);

}

#endif
