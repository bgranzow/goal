#include <cassert>

#include <apf.h>
#include <apfCavityOp.h>
#include <apfMesh.h>
#include <PCU.h>

#include "goal_control.hpp"
#include "goal_size_field.hpp"

namespace goal {

struct Specification {
  apf::Mesh* mesh;
  apf::Field* vtx_error;
  apf::Field* error;
  int polynomial_order;
  size_t target_number;
  double alpha;
  double beta;
  double size_factor;
  apf::Field* elem_size;
  apf::Field* vtx_size;
};

static void setup_specification(
    Specification* s, apf::Field* err, size_t t, int p) {
  s->mesh = apf::getMesh(err);
  s->vtx_error = err;
  s->error = 0;
  s->polynomial_order = p;
  s->target_number = t;
  s->alpha = 0.25;
  s->beta = 2.0;
  s->size_factor = 0.0;
  s->elem_size = 0;
  s->vtx_size = 0;
}

static void interpolate_error(Specification* s) {
  auto m = s->mesh;
  s->error = apf::createStepField(s->mesh, "elem_e", apf::SCALAR);
  apf::MeshEntity* elem;
  apf::Vector3 p;
  auto it = m->begin(m->getDimension());
  while ((elem = s->mesh->iterate(it))) {
    auto me = apf::createMeshElement(m, elem);
    auto fe = apf::createElement(s->vtx_error, me);
    apf::getIntPoint(me, 1, 0, p);
    auto v = apf::getScalar(fe, p);
    apf::setScalar(s->error, elem, 0, v);
    apf::destroyElement(fe);
    apf::destroyMeshElement(me);
  }
}

static double sum_contributions(Specification* s) {
  double r = 0.0;
  auto d = s->mesh->getDimension();
  auto p = s->polynomial_order;
  apf::MeshEntity* elem;
  auto it = s->mesh->begin(d);
  while ((elem = s->mesh->iterate(it))) {
    auto v = fabs(apf::getScalar(s->error, elem, 0));
    r += pow(v, ((2.0 * d) / (2.0 * p + d)));
  }
  s->mesh->end(it);
  PCU_Add_Doubles(&r, 1);
  return r;
}

static void compute_size_factor(Specification* s) {
  auto d = s->mesh->getDimension();
  double G = sum_contributions(s);
  double N = s->target_number;
  s->size_factor = pow((G / N), (1.0 / d));
}

static double get_current_size(apf::Mesh* m, apf::MeshEntity* e) {
  double h = 0.0;
  apf::Downward edges;
  int ne = m->getDownward(e, 1, edges);
  for (int i = 0; i < ne; ++i) h = std::max(h, apf::measure(m, edges[i]));
  return h;
}

static double get_new_size(Specification* s, apf::MeshEntity* e) {
  auto p = s->polynomial_order;
  auto d = s->mesh->getDimension();
  auto h = get_current_size(s->mesh, e);
  auto theta_e = fabs(apf::getScalar(s->error, e, 0));
  auto r = pow(theta_e, ((-2.0) / (2.0 * p + d)));
  auto h_new = s->size_factor * r * h;
  if (h_new < s->alpha * h) h_new = s->alpha * h;
  if (h_new > s->beta * h) h_new = s->beta * h;
  return h_new;
}

static void get_elem_size(Specification* s) {
  auto e_size = apf::createStepField(s->mesh, "esize", apf::SCALAR);
  auto d = s->mesh->getDimension();
  apf::MeshEntity* elem;
  auto it = s->mesh->begin(d);
  while ((elem = s->mesh->iterate(it))) {
    auto h = get_new_size(s, elem);
    apf::setScalar(e_size, elem, 0, h);
  }
  s->mesh->end(it);
  s->elem_size = e_size;
}

static void avg_to_vtx(apf::Field* ef, apf::Field* vf, apf::MeshEntity* ent) {
  auto m = apf::getMesh(ef);
  apf::Adjacent elems;
  m->getAdjacent(ent, m->getDimension(), elems);
  double s = 0.0;
  for (size_t i = 0; i < elems.getSize(); ++i)
    s += apf::getScalar(ef, elems[i], 0);
  s /= elems.getSize();
  apf::setScalar(vf, ent, 0, s);
}

class AverageOp : public apf::CavityOp {
 public:
  AverageOp(Specification* s) : apf::CavityOp(s->mesh), specs(s), entity(0) {}
  virtual Outcome setEntity(apf::MeshEntity* e) {
    entity = e;
    if (apf::hasEntity(specs->vtx_size, entity)) return SKIP;
    if (!requestLocality(&entity, 1)) return REQUEST;
    return OK;
  }
  virtual void apply() {
    avg_to_vtx(specs->elem_size, specs->vtx_size, entity);
  }
  Specification* specs;
  apf::MeshEntity* entity;
};

static void average_size_field(Specification* s) {
  s->vtx_size = apf::createLagrangeField(s->mesh, "size", apf::SCALAR, 1);
  AverageOp op(s);
  op.applyToDimension(0);
}

static void create_size_field(Specification* s) {
  interpolate_error(s);
  compute_size_factor(s);
  get_elem_size(s);
  average_size_field(s);
  apf::destroyField(s->elem_size);
  apf::destroyField(s->error);
}

apf::Field* get_iso_target_size(apf::Field* e, std::size_t t, int p) {
  auto t0 = time();
  assert(t > 0);
  Specification s;
  setup_specification(&s, e, t, p);
  create_size_field(&s);
  auto t1 = time();
  print("isotropic target size field computed in %f seconds", t1 - t0);
  return s.vtx_size;
}

} /* namespace goal */
