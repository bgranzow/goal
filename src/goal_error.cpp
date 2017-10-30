#include <apf.h>
#include <apfMesh.h>
#include <PCU.h>

namespace goal {

apf::Field* compute_error(apf::Field* ue, apf::Field* pe) {
  double p_err = 0.0;
  apf::Vector3 u_err(0,0,0);
  apf::Vector3 xi(0,0,0);
  auto m = apf::getMesh(ue);
  int num_dims = m->getDimension();
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = m->begin(num_dims);
  auto err = apf::createStepField(m, "err", apf::SCALAR);
  while ((elem = m->iterate(elems))) {
    auto me = apf::createMeshElement(m, elem);
    auto u_elem = apf::createElement(ue, me);
    auto p_elem = apf::createElement(pe, me);
    apf::getIntPoint(me, 1, 0, xi);
    p_err = apf::getScalar(p_elem, xi);
    apf::getVector(u_elem, xi, u_err);
    double total = 0.0;
    for (int d = 0; d < num_dims; ++d)
      total += u_err[d];
    total = std::abs(total);
    total += std::abs(p_err);
    apf::setScalar(err, elem, 0, total);
    apf::destroyElement(u_elem);
    apf::destroyElement(p_elem);
    apf::destroyMeshElement(me);
  }
  m->end(elems);
  return err;
}

double sum_contribs(apf::Field* e) {
  auto m = apf::getMesh(e);
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = m->begin(m->getDimension());
  double sum = 0.0;
  while ((elem = m->iterate(elems)))
    sum += apf::getScalar(e, elem, 0);
  m->end(elems);
  PCU_Add_Doubles(&sum, 1);
  return sum;
}

}
