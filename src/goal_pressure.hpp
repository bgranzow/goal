#ifndef goal_pressure_hpp
#define goal_pressure_hpp

#include <apf.h>
#include "goal_integrator.hpp"
#include "goal_scalar_types.hpp"

namespace goal {

class Disc;
template <typename T> class Pressure;

template <>
class Pressure<ST> : public Integrator {
  public:
    Pressure(apf::Field* base, const int mode);
    ~Pressure();
    int get_num_dims() { return num_dims; }
    int get_num_nodes() { return num_nodes; }
    ST& val();
    ST& grad(const int i);
    ST& nodal(const int n);
    ST& resid(const int n);
    void pre_process(SolInfo* s);
    void gather(apf::MeshElement* me);
    void at_point(apf::Vector3 const& p, double, double);
    void scatter(SolInfo* s);
    void post_process(SolInfo*);
  private:
    std::function<void(Pressure<ST>*, SolInfo*)> op;
    void scatter_none(SolInfo* s);
    void scatter_primal(SolInfo* s);
    Disc* disc;
    apf::Field* field;
    apf::FieldShape* shape;
    apf::Element* elem;
    apf::NewArray<ST> BF;
    apf::NewArray<apf::Vector3> GBF;
    apf::NewArray<ST> node;
    ST value;
    apf::Vector3 gradient;
    std::vector<ST> residual;
    int num_dims;
    int num_nodes;
};

template <>
class Pressure<FADT> : public Integrator {
  public:
    Pressure(apf::Field* base, const int mode);
    ~Pressure();
    int get_num_dims() { return num_dims; }
    int get_num_nodes() { return num_nodes; }
    FADT& val();
    FADT& grad(const int i);
    FADT& nodal(const int n);
    FADT& resid(const int n);
    void pre_process(SolInfo* s);
    void gather(apf::MeshElement* me);
    void at_point(apf::Vector3 const& p, double, double);
    void scatter(SolInfo* s);
    void post_process(SolInfo*);
  private:
    std::function<void(Pressure<FADT>*, SolInfo*)> op;
    void scatter_none(SolInfo* s);
    void scatter_primal(SolInfo* s);
    void scatter_adjoint(SolInfo* s);
    Disc* disc;
    apf::Field* field;
    apf::FieldShape* shape;
    apf::Element* elem;
    apf::NewArray<ST> BF;
    apf::NewArray<apf::Vector3> GBF;
    apf::NewArray<ST> node_st;
    std::vector<FADT> node_fadt;
    FADT value;
    std::vector<FADT> gradient;
    std::vector<FADT> residual;
    int num_dims;
    int num_nodes;
    int num_dofs;
};

}

#endif
