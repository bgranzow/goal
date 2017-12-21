#ifndef goal_point_wise_hpp
#define goal_point_wise_hpp

#include <Teuchos_ParameterList.hpp>
#include "goal_qoi.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

template <typename T>
class PointWise : public QoI<T> {
  public:
    PointWise(ParameterList const& p);
    void pre_process(SolInfo* s);
    void post_process(SolInfo* s);
  private:
    ParameterList params;
};

}

#endif
