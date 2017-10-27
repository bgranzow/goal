#ifndef goal_tbcs_hpp
#define goal_tbcs_hpp

#include <Teuchos_ParameterList.hpp>

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class SolInfo;
class Integrator;

void set_tbcs(
    ParameterList const& p,
    RCP<Integrator> w,
    SolInfo* s,
    double t);

}

#endif
