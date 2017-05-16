#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_ev_dual.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename EVALT, typename TRAITS>
Dual<EVALT, TRAITS>::Dual(
    std::vector<Field*> const& z_,
    std::vector<Field*> const& z_fine_,
    int type)
    : z(z_),
      z_fine(z_fine_) {
  GOAL_DEBUG_ASSERT(z.size() == z_fine.size());
  num_fields = z.size();
  num_vtx = z_fine[0]->get_num_vtx(type);
  num_ips = z_fine[0]->get_num_ips(type);
  num_dims = z_fine[0]->get_num_dims();
  w.resize(num_fields);
  gw.resize(num_fields);
  for (int i = 0; i < num_fields; ++i) {
    auto name = z_fine[i]->name();
    auto gname = z_fine[i]->g_name();
    auto dl = z_fine[i]->PU_dl(type);
    auto gdl = z_fine[i]->g_PU_dl(type);
    w[i] = PHX::MDField<ScalarT, Ent, Node>(name, dl);
    gw[i] = PHX::MDField<ScalarT, Ent, Node, Dim>(gname, gdl);
    this->addEvaluatedField(w[i]);
    this->addEvaluatedField(gw[i]);
  }
  this->setName("Dual");
}

PHX_POST_REGISTRATION_SETUP(Dual, data, fm) {
  for (int i = 0; i < num_fields; ++i) {
    this->utils.setFieldData(w[i], fm);
    this->utils.setFieldData(gw[i], fm);
  }
}

PHX_EVALUATE_FIELDS(Dual, workset) {
  apf::Vector3 p;
  apf::Vector3 grad_z_val;
  apf::Vector3 grad_z_fine_val;
  apf::Vector3 grad_z_fine_min_z_val;
  apf::NewArray<double> PU;
  apf::NewArray<apf::Vector3> GPU;
  auto m = z_fine[0]->get_apf_mesh();
  auto q = z_fine[0]->get_q_degree();
  auto s = apf::getLagrange(1);
  for (int elem = 0; elem < workset.size; ++elem) {
    auto e = workset.entities[elem];
    auto me = apf::createMeshElement(m, e);
    for (int field = 0; field < num_fields; ++field) {
      auto apf_z = z[field]->get_apf_field();
      auto apf_z_fine = z_fine[field]->get_apf_field();
      auto z_elem = apf::createElement(apf_z, me);
      auto z_fine_elem = apf::createElement(apf_z_fine, me);
      for (int ip = 0; ip < num_ips; ++ip) {
        apf::getIntPoint(me, q, ip, p);
        apf::getBF(s, me, p, PU);
        apf::getGradBF(s, me, p, GPU);
        auto z_val = apf::getScalar(z_elem, p);
        auto z_fine_val = apf::getScalar(z_fine_elem, p);
        auto z_fine_min_z_val = z_fine_val - z_val;
        apf::getGrad(z_elem, p, grad_z_val);
        apf::getGrad(z_fine_elem, p, grad_z_fine_val);
        grad_z_fine_min_z_val = grad_z_fine_val - grad_z_val;
        for (int vtx = 0; vtx < num_vtx; ++vtx) {
          w(elem, vtx, ip) = z_fine_min_z_val * PU[vtx];
          for (int dim = 0; dim < num_dims; ++dim)
            gw(elem, vtx, ip, dim) =
              grad_z_fine_min_z_val[dim] * PU[vtx] +
              z_fine_min_z_val * GPU[vtx][dim];
        }
      }
      apf::destroyElement(z_elem);
      apf::destroyElement(z_fine_elem);
    }
    apf::destroyMeshElement(me);
  }
}

} // end namespace goal
