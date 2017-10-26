#ifndef goal_mixed_hpp
#define goal_mixed_hpp

#include <Teuchos_RCP.hpp>
#include "goal_integrator.hpp"

namespace goal {

using Teuchos::RCP;

class States;
template <typename T> class Model;
template <typename T> class Pressure;

template <typename T>
class Mixed : public Integrator {
  public:
    Mixed(
        RCP<Integrator> p,
        RCP<Model<T>> cm,
        States* s,
        const bool save);
    void in_elem(apf::MeshElement* me);
    void at_point(apf::Vector3 const&, double, double);
    void out_elem();
  private:
    RCP<Pressure<T>> p;
    RCP<Model<T>> model;
    States* states;
    bool save_state;
    apf::MeshEntity* elem;
    int num_dims;
};

}

#endif
