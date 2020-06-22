#ifndef goal_assembly_hpp
#define goal_assembly_hpp

#include <vector>
#include <Teuchos_RCP.hpp>

namespace goal {

class SolInfo;
class Integrator;

using Teuchos::RCP;
using Evaluators = std::vector<RCP<Integrator>>;

RCP<Integrator> find_evaluator(std::string const& n, Evaluators const& E);
void set_time(Evaluators& E, double t_now, double t_old);
void assemble(Evaluators const& E, SolInfo* s);

}

#endif
