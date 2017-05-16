#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_ev_basis.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename EVALT, typename TRAITS>
Basis<EVALT, TRAITS>::Basis(Field* f, int type)
    : field(f),
      num_nodes(f->get_num_nodes(type)),
      num_ips(f->get_num_ips(type)),
      num_dims(f->get_num_dims()),
      wdv(f->wdv_name(), f->ip0_dl(type)),
      bf(f->basis_name(), f->w_dl(type)),
      gbf(f->g_basis_name(), f->g_w_dl(type)) {
  this->addEvaluatedField(wdv);
  this->addEvaluatedField(bf);
  this->addEvaluatedField(gbf);
  this->setName(f->basis_name());
}

PHX_POST_REGISTRATION_SETUP(Basis, data, fm) {
  this->utils.setFieldData(wdv, fm);
  this->utils.setFieldData(bf, fm);
  this->utils.setFieldData(gbf, fm);
  (void)data;
}

PHX_EVALUATE_FIELDS(Basis, workset) {
  apf::Vector3 p;
  apf::NewArray<double> BF;
  apf::NewArray<apf::Vector3> GBF;
  auto q = field->get_q_degree();
  auto s = field->get_apf_basis();
  auto m = field->get_apf_mesh();
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto me = apf::createMeshElement(m, e);
    for (int ip = 0; ip < num_ips; ++ip) {
      apf::getIntPoint(me, q, ip, p);
      auto w = apf::getIntWeight(me, q, ip);
      wdv(elem, ip) = w * apf::getDV(me, p);
      apf::getBF(s, me, p, BF);
      apf::getGradBF(s, me, p, GBF);
      for (int node = 0; node < num_nodes; ++node) {
        bf(elem, node, ip) = BF[node];
        for (int dim = 0; dim < num_dims; ++dim)
          gbf(elem, node, ip, dim) = GBF[node][dim];
      }
    }
    apf::destroyMeshElement(me);
  }
}

template class Basis<goal::Traits::Residual, goal::Traits>;
template class Basis<goal::Traits::Jacobian, goal::Traits>;

} // end namespace goal
