#include "goal_control.hpp"
#include "goal_nested.hpp"

#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>

namespace goal {

Nested::Nested(Disc* d, int m) {
  double t0 = time();
  mode = m;
  am_base = false;
  sets = d->get_model_sets();
  base_mesh = d->get_apf_mesh();
  number_elems();
  copy_mesh();
  initialize();
  double t1 = time();
  print(" > nested mesh built in %f seconds", t1 - t0);
}

Nested::~Nested() {
  (void)old_vtx_tag;
  (void)new_vtx_tag;
}

void Nested::number_elems() {
  int i = 0;
  int dim = base_mesh->getDimension();
  base_elems.resize(base_mesh->count(dim));
  auto vt = apf::SCALAR;
  auto s = apf::getIPFitShape(dim, 1);
  auto f = apf::createField(base_mesh, "elems", vt, s);
  apf::MeshEntity* elem;
  auto it = base_mesh->begin(dim);
  while ((elem = base_mesh->iterate(it))) {
    apf::setScalar(f, elem, 0, i++);
    base_elems[i] = elem;
  }
  base_mesh->end(it);
}

void Nested::copy_mesh() {
  auto model = base_mesh->getModel();
  mesh = apf::createMdsMesh(model, base_mesh);
  apf::disownMdsModel(mesh);
}

Nested* create_nested(Disc* d, int m) {
  return new Nested(d, m);
}

void destroy_nested(Nested* n) {
  delete n;
}

}
