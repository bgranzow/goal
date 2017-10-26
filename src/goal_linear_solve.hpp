#ifndef goal_linear_solve_hpp
#define goal_linear_solve_hpp

#include "goal_data_types.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Disc;

void solve(
    ParameterList const& p,
    RCP<MatrixT> A,
    RCP<VectorT> x,
    RCP<VectorT> b,
    Disc* d);

}

#endif
