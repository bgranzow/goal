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
    void localize(double t_now, double t_old);
    ParameterList params;
    Primal* primal;
    Disc* base_disc;
    Nested* nested_disc;
    Mechanics* mech;
    SolInfo* sol_info;
    Evaluators adjoint;
    Evaluators error;
    apf::Field* z_disp;
    apf::Field* z_press;
    apf::Field* z_disp_diff;
    apf::Field* z_press_diff;
    apf::Field* e_disp;
    apf::Field* e_press;
};

NestedAdjoint* create_nested_adjoint(ParameterList const& p, Primal* pr);
void destroy_nested_adjoint(NestedAdjoint* a);

}

#endif
