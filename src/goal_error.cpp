#include <apf.h>
#include <apfMesh.h>
#include <PCU.h>

namespace goal {

apf::Field* compute_error(apf::Field* ue) {
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
    apf::getIntPoint(me, 1, 0, xi);
    apf::getVector(u_elem, xi, u_err);
    double total = 0.0;
    for (int d = 0; d < num_dims; ++d)
      total += std::abs(u_err[d]);
    apf::setScalar(err, elem, 0, total);
    apf::destroyElement(u_elem);
    apf::destroyMeshElement(me);
  }
  m->end(elems);
  return err;
}

double sum_contribs(apf::Field* ue) {
  double sum = 0.0;
  apf::Vector3 u(0,0,0);
  apf::MeshEntity* vtx;
  auto m = apf::getMesh(ue);
  auto it = m->begin(0);
  auto dim = m->getDimension();
  while ((vtx = m->iterate(it))) {
    apf::getVector(ue, vtx, 0, u);
    for (int d = 0; d < dim; ++d)
      sum += std::abs(u[d]);
  }
  m->end(it);
  PCU_Add_Doubles(&sum, 1);
  return sum;
}

}
