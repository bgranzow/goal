#ifndef goal_mechanics_hpp
#define goal_mechanics_hpp

#include <Teuchos_ParameterList.hpp>

namespace apf {
class Field;
}

namespace goal {

class Disc;
class Integrator;
class SolInfo;
class States;

using Teuchos::RCP;
using Teuchos::ParameterList;
using Evaluators = std::vector<RCP<Integrator>>;

class Mechanics {
  public:
    Mechanics(ParameterList const& p, Disc* d);
    ~Mechanics();
    Disc* get_disc() { return disc; }
    States* get_states() { return states; }
    apf::Field* get_displacement() { return displacement; }
    template <typename T>
    void build_resid(Evaluators& E, bool save_states);
    template <typename T>
    void build_functional(ParameterList const& params, Evaluators& E);
    void build_error(Evaluators& E);
  private:
    void make_displacement();
    void make_states();
    ParameterList params;
    Disc* disc;
    States* states;
    std::string model;
    apf::Field* displacement;
};

Mechanics* create_mechanics(ParameterList const& p, Disc* d);
void destroy_mechanics(Mechanics* m);

}

#endif
