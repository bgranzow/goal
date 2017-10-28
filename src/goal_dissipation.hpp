#ifndef goal_dissipation_hpp
#define goal_dissipation_hpp

#include <Teuchos_ParameterList.hpp>

namespace apf {
class Field;
}

namespace goal {

class Disc;
class Integrator;
class Mechanics;
class Primal;
class SolInfo;

using Teuchos::RCP;
using Teuchos::ParameterList;
using Evaluators = std::vector<RCP<Integrator>>;

class Dissipation {
  public:
    Dissipation(ParameterList const& p, Primal* pr);
    ~Dissipation();
    void run(double t_now, double t_old);;
  private:
    void print_banner(double t_now);
    void build_adjoint();
    void destroy_adjoint();
    void compute_adjoint(double t_now, double t_old);
    void solve_adjoint(double t_now, double t_old);
    ParameterList params;
    Primal* primal;
    Disc* disc;
    Mechanics* mech;
    SolInfo* sol_info;
    Evaluators adjoint;
    apf::Field* z_displacement;
    apf::Field* z_pressure;
};

Dissipation* create_dissipation(ParameterList const& p, Primal* pr);
void destroy_dissipation(Dissipation* d);

}

#endif
