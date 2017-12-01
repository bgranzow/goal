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
    total += p_err;
    total = std::abs(total);
    apf::setScalar(err, elem, 0, total);
    apf::destroyElement(u_elem);
    apf::destroyElement(p_elem);
    apf::destroyMeshElement(me);
  }
  m->end(elems);
  return err;
}

double sum_contribs(apf::Field* ue, apf::Field* pe) {
  double sum = 0.0;
  apf::Vector3 u(0,0,0);
  apf::MeshEntity* vtx;
  auto m = apf::getMesh(ue);
  auto it = m->begin(0);
  auto dim = m->getDimension();
  while ((vtx = m->iterate(it))) {
    apf::getVector(ue, vtx, 0, u);
    auto p = apf::getScalar(pe, vtx, 0);
    double tmp = 0.0;
    for (int d = 0; d < dim; ++d)
      tmp += u[d];
    tmp += p;
    sum += std::abs(tmp);
  }
  m->end(it);
  PCU_Add_Doubles(&sum, 1);
  return sum;
}

}
