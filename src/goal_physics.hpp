#ifndef goal_physics_hpp
#define goal_physics_hpp

#include <Teuchos_ParameterList.hpp>

namespace apf {
class Field;
}

namespace goal {

class Disc;
class Integrator;
class SolInfo;

using Teuchos::RCP;
using Teuchos::ParameterList;
using Evaluators = std::vector<RCP<Integrator>>;

class Physics {
  public:
    Physics(ParameterList const& p, Disc* d);
    ~Physics();
    Disc* get_disc() { return disc; }
    apf::Field* get_soln() { return soln; }
    template <typename T>
    void build_resid(Evaluators& E);
    template <typename T>
    void build_functional(ParameterList const& params, Evaluators& E);
    void build_error(Evaluators& E);
  private:
    void make_soln();
    void make_states();
    ParameterList params;
    Disc* disc;
    std::string model;
    apf::Field* soln;
};

Physics* create_physics(ParameterList const& p, Disc* d);
void destroy_physics(Physics* m);

}

#endif
