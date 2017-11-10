#ifndef goal_ibcs_hpp
#define goal_ibcs_hpp

#include <Teuchos_ParameterList.hpp>

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class SolInfo;
class Integrator;

void set_ibcs(
    ParameterList const& p,
    RCP<Integrator> w,
    SolInfo* s,
    double t);

}

#endif
