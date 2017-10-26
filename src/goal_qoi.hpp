#ifndef goal_qoi_hpp
#define goal_qoi_hpp

#include "goal_integrator.hpp"
#include "goal_scalar_types.hpp"

namespace goal {

class Disc;
template <typename T> class QoI;

template <>
class QoI<ST> : public Integrator {
  public:
    QoI();
    virtual ~QoI();
    ST const& get_qoi_value() const { return qoi_value; }
    ST const& get_elem_value() const { return elem_value; }
    virtual void set_time(const double, const double) {}
    virtual void set_elem_set(const int) {}
    virtual void pre_process(SolInfo*);
    virtual void gather(apf::MeshElement* me);
    virtual void in_elem(apf::MeshElement*) {}
    virtual void at_point(apf::Vector3 const&, double, double) {}
    virtual void out_elem() {}
    virtual void scatter(SolInfo* s);
    virtual void post_process(SolInfo*);
  protected:
    apf::MeshElement* elem;
    Disc* disc;
    ST qoi_value;
    ST elem_value;
};

template <>
class QoI<FADT> : public Integrator {
  public:
    QoI();
    virtual ~QoI();
    ST const& get_qoi_value() const { return qoi_value; }
    FADT const& get_elem_value() const { return elem_value; }
    virtual void set_time(const double, const double) {}
    virtual void set_elem_set(const int) {}
    virtual void pre_process(SolInfo* s);
    virtual void gather(apf::MeshElement* me);
    virtual void in_elem(apf::MeshElement*) {}
    virtual void at_point(apf::Vector3 const&, double, double) {}
    virtual void out_elem() {}
    virtual void scatter(SolInfo* s);
    virtual void post_process(SolInfo* s);
  protected:
    apf::MeshElement* elem;
    Disc* disc;
    ST qoi_value;
    FADT elem_value;
};

}

#endif
