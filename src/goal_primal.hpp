#ifndef goal_primal_hpp
#define goal_primal_hpp

#include <Teuchos_ParameterList.hpp>

namespace goal {

class Integrator;
class Mechanics;
class SolInfo;

using Teuchos::RCP;
using Teuchos::ParameterList;
using Evaluators = std::vector<RCP<Integrator>>;

class Primal {
  public:
    Primal(ParameterList const& p, Mechanics* m);
    ~Primal();
    Mechanics* get_mech() { return mech; }
    SolInfo* get_sol_info() { return sol_info; }
    void build_data();
    void destroy_data();
    void solve(const double t_now, const double t_old);
  private:
    void print_banner(const double t_now);
    void compute_resid(const double t_now, const double t_old);
    void compute_jacob(const double t_now, const double t_old);
    ParameterList params;
    Mechanics* mech;
    SolInfo* sol_info;
    Evaluators residual;
    Evaluators jacobian;
};

Primal* create_primal(ParameterList const& p, Mechanics* m);
void destroy_primal(Primal* p);

}

#endif
