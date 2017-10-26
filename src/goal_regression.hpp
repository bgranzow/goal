#ifndef goal_regression_hpp
#define goal_regression_hpp

namespace Teuchos {
class ParameterList;
}

namespace goal {

using Teuchos::ParameterList;

class Functional;

void check_J_regression(ParameterList const& p, Functional* f);

}

#endif
