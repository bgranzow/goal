#ifndef goal_displacement_hpp
#define goal_displacement_hpp

#include <apf.h>
#include "goal_integrator.hpp"
#include "goal_scalar_types.hpp"

namespace goal {

class Disc;
template <typename T> class Displacement;

template <>
class Displacement<ST> : public Integrator {
  public:
    Displacement(apf::Field* base, int mode);
    ~Displacement();
    int get_num_dims() { return num_dims; }
    int get_num_nodes() { return num_nodes; }
    ST& val(int i);
    ST& grad(int i, int j);
    ST& nodal(int n, int i);
    ST& resid(int n, int i);
    void pre_process(SolInfo* s);
    void gather(apf::MeshElement* me);
    void at_point(apf::Vector3 const& p, double, double);
    void scatter(SolInfo* s);
    void post_process(SolInfo*);
  private:
    std::function<void(Displacement<ST>*, SolInfo*)> op;
    void scatter_none(SolInfo* s);
    void scatter_primal(SolInfo* s);
    Disc* disc;
    apf::Field* field;
    apf::FieldShape* shape;
    apf::Element* elem;
    apf::NewArray<ST> BF;
    apf::NewArray<apf::Vector3> GBF;
    apf::NewArray<apf::Vector3> node;
    apf::Vector3 value;
    apf::Matrix3x3 gradient;
    std::vector<ST> residual;
    int num_dims;
    int num_nodes;
};

template <>
class Displacement<FADT> : public Integrator {
  public:
    Displacement(apf::Field* base, int mode);
    ~Displacement();
    int get_num_dims() { return num_dims; }
    int get_num_nodes() { return num_nodes; }
    FADT& val(int i);
    FADT& grad(int i, int j);
    FADT& nodal(int n, int i);
    FADT& resid(int n, int i);
    void pre_process(SolInfo* s);
    void gather(apf::MeshElement* me);
    void at_point(apf::Vector3 const& p, double, double);
    void scatter(SolInfo* s);
    void post_process(SolInfo*);
  private:
    std::function<void(Displacement<FADT>*, SolInfo*)> op;
    void scatter_none(SolInfo* s);
    void scatter_primal(SolInfo* s);
    void scatter_adjoint(SolInfo* s);
    Disc* disc;
    apf::Field* field;
    apf::FieldShape* shape;
    apf::Element* elem;
    apf::NewArray<ST> BF;
    apf::NewArray<apf::Vector3> GBF;
    apf::NewArray<apf::Vector3> node_st;
    std::vector<FADT> node_fadt;
    std::vector<FADT> value;
    std::vector<FADT> gradient;
    std::vector<FADT> residual;
    int num_dims;
    int num_nodes;
    int num_dofs;
};

}

#endif
