#include <apf.h>
#include <apfMesh.h>
#include <PCU.h>

namespace goal {

apf::Field* compute_error(apf::Field* e) {
  apf::Vector3 xi(0,0,0);
  auto m = apf::getMesh(e);
  int num_dims = m->getDimension();
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = m->begin(num_dims);
  auto err = apf::createStepField(m, "err", apf::SCALAR);
  while ((elem = m->iterate(elems))) {
    auto me = apf::createMeshElement(m, elem);
    auto e_elem = apf::createElement(e, me);
    apf::getIntPoint(me, 1, 0, xi);
    auto elem_err = std::abs(apf::getScalar(e_elem, xi));
    apf::setScalar(err, elem, 0, elem_err);
    apf::destroyElement(e_elem);
    apf::destroyMeshElement(me);
  }
  m->end(elems);
  return err;
}

double sum_contribs(apf::Field* e) {
  double sum = 0.0;
  apf::MeshEntity* vtx;
  auto m = apf::getMesh(e);
  auto it = m->begin(0);
  while ((vtx = m->iterate(it)))
    sum += std::abs(apf::getScalar(e, vtx, 0));
  m->end(it);
  PCU_Add_Doubles(&sum, 1);
  return sum;
}

}
