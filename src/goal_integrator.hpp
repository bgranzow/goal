#ifndef goal_integrator_hpp
#define goal_integrator_hpp

#include <string>

namespace apf {
class MeshEntity;
class Vector3;
class VectorElement;
using MeshElement = VectorElement;
}

namespace goal {

class SolInfo;

class Integrator {
  public:
    Integrator();
    virtual ~Integrator();
    std::string const& get_name() { return name; }
    virtual void set_time(const double, const double) {}
    virtual void pre_process(SolInfo*) {}
    virtual void set_elem_set(const int) {}
    virtual void gather(apf::MeshElement*) {}
    virtual void in_elem(apf::MeshElement*) {}
    virtual void at_point(apf::Vector3 const&, double, double) {}
    virtual void out_elem() {}
    virtual void scatter(SolInfo*) {}
    virtual void post_process(SolInfo*) {}
  protected:
    std::string name;
};

}

#endif
