#ifndef goal_ks_vm_hpp
#define goal_ks_vm_hpp

#include <Teuchos_ParameterList.hpp>
#include "goal_qoi.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

class Integrator;
class SolInfo;
template <typename T> class Model;

template <typename T>
class KSVM : public QoI<T> {
  public:
    KSVM(ParameterList const& p, RCP<Integrator> model);
    void pre_process(SolInfo* s);
    void at_point(apf::Vector3 const&, double w, double dv);
    void post_process(SolInfo* s);
  private:
    ParameterList params;
    RCP<Model<T>> model;
    double max;
    double rho;
    double scale;
    int num_dims;
};

}

#endif
