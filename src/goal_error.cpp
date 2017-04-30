#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <PCU.h>

#include "goal_error.hpp"
#include "goal_field.hpp"

namespace goal {

static void validate_fields(std::vector<RCP<Field> > const& e) {
  for (size_t i = 0; i < e.size(); ++i) {
    assert(e[i]->get_basis_type() == LAGRANGE);
    assert(e[i]->get_p_order() == 1);
  }
}

class Iterator {
 public:
  Iterator(std::vector<RCP<Field> > const& error) {
    e = error;
    m = e[0]->get_apf_mesh();
    vtx = 0;
    vals[0] = 0.0; vals[1] = 0.0; vals[2] = 0.0;
  }
  void apply() {
    apf::MeshIterator* vertices = m->begin(0);
    while ((vtx = m->iterate(vertices))) {
      if (! m->isOwned(vtx)) continue;
      for (size_t f = 0; f < e.size(); ++f) {
        auto apf_f = e[f]->get_apf_field();
        int nc = apf::countComponents(apf_f);
        apf::getComponents(apf_f, vtx, 0, vals);
        pre_component_op();
        for (int c = 0; c < nc; ++c)
          component_op(c);
        post_component_op();
      }
    }
    m->end(vertices);
    synchronize();
  }
  virtual void pre_component_op() = 0;
  virtual void component_op(int c) = 0;
  virtual void post_component_op() = 0;
  virtual void synchronize() = 0;

 protected:
  apf::Mesh* m;
  apf::MeshEntity* vtx;
  std::vector<RCP<Field> > e;
  double vals[3];
};

class Adder : public Iterator {
 public:
  Adder(std::vector<RCP<Field> > const& e)
    : Iterator(e),
      sum(0.0) {}
  void pre_component_op() {}
  void component_op(int c) { sum += vals[c]; }
  void post_component_op() {}
  void synchronize() { PCU_Add_Doubles(&sum, 1); }
  double sum;
};

class Bounder : public Iterator {
 public:
  Bounder(std::vector<RCP<Field> > const& e)
    : Iterator(e),
      bound(0.0) {}
  void pre_component_op() { tmp = 0.0; }
  void component_op(int c) { tmp += vals[c]; }
  void post_component_op() { bound += std::abs(tmp); }
  void synchronize() { PCU_Add_Doubles(&bound, 1); }
  double bound;

 private:
  double tmp;
};

class Indicator : public Iterator {
 public:
  Indicator(std::vector<RCP<Field> > const& e)
    : Iterator(e) {
    ind = apf::createLagrangeField(m, "ind", apf::SCALAR, 1);
  }
  void pre_component_op() { tmp = 0.0; }
  void component_op(int c) { tmp += vals[c]; }
  void post_component_op() { apf::setScalar(ind, vtx, 0, std::abs(tmp)); }
  void synchronize() { apf::synchronize(ind); }
  apf::Field* ind;

 private:
  double tmp;
};

double sum_contributions(std::vector<RCP<Field> > const& e) {
  validate_fields(e);
  Adder adder(e);
  adder.apply();
  return adder.sum;
}

double approx_upper_bound(std::vector<RCP<Field> > const& e) {
  validate_fields(e);
  Bounder bounder(e);
  bounder.apply();
  return bounder.bound;
}

apf::Field* compute_indicators(std::vector<RCP<Field> > const& e) {
  validate_fields(e);
  Indicator indicator(e);
  indicator.apply();
  return indicator.ind;
}

} /* namespace goal */
