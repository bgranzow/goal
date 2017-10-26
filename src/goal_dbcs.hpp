#ifndef goal_dbcs_hpp
#define goal_dbcs_hpp

namespace Teuchos {
class ParameterList;
}

namespace goal {

using Teuchos::ParameterList;

class SolInfo;

void set_resid_dbcs(ParameterList const& p, SolInfo* s, const double t);
void set_jac_dbcs(ParameterList const& p, SolInfo* s, const double t);

}

#endif
