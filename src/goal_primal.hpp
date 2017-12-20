#ifndef goal_primal_hpp
#define goal_primal_hpp

#include <Teuchos_ParameterList.hpp>

namespace goal {

class Integrator;
class Physics;
class SolInfo;

using Teuchos::RCP;
using Teuchos::ParameterList;
using Evaluators = std::vector<RCP<Integrator>>;

class Primal {
  public:
    Primal(ParameterList const& p, Physics* phy);
    ~Primal();
    Physics* get_physics() { return physics; }
    SolInfo* get_sol_info() { return sol_info; }
    void build_data();
    void destroy_data();
    void solve(double t_now, double t_old);
  private:
    void print_banner(double t_now);
    void compute_resid(double t_now, double t_old);
    void compute_jacob(double t_now, double t_old);
    ParameterList params;
    Physics* physics;
    SolInfo* sol_info;
    Evaluators residual;
    Evaluators jacobian;
};

Primal* create_primal(ParameterList const& p, Physics* phy);
void destroy_primal(Primal* p);

}

#endif
