#include "goal_ev_interpolate.hpp"
#include "goal_field.hpp"
#include "goal_traits.hpp"
#include "goal_workset.hpp"

namespace goal {

template <typename EVALT, typename TRAITS>
Interpolate<EVALT, TRAITS>::Interpolate(
    std::vector<Field*> const& f, int type)
    : bf(f[0]->basis_name(), f[0]->w_dl(type)),
      gbf(f[0]->g_basis_name(), f[0]->g_w_dl(type)) {
  num_fields = f.size();
  num_nodes = f[0]->get_num_nodes(type);
  num_ips = f[0]->get_num_ips(type);
  num_dims = f[0]->get_num_dims();
  nodal.resize(num_fields);
  u.resize(num_fields);
  gu.resize(num_fields);
  for (int i = 0; i < num_fields; ++i) {
    auto name = f[i]->name();
    auto gname = f[i]->g_name();
    auto dl = f[i]->dl(type);
    auto ip_dl = f[i]->ip_dl(type);
    auto g_ip_dl = f[i]->g_ip_dl(type);
    nodal[i] = PHX::MDField<ScalarT, Ent, Node>(name, dl);
    u[i] = PHX::MDField<ScalarT, Ent, Node>(name, ip_dl);
    gu[i] = PHX::MDField<ScalarT, Ent, Node, Dim>(gname, g_ip_dl);
    this->addDependentField(nodal[i]);
    this->addEvaluatedField(u[i]);
    this->addEvaluatedField(gu[i]);
  }
  this->addDependentField(bf);
  this->addDependentField(gbf);
  this->setName("Interpolate");
}

PHX_POST_REGISTRATION_SETUP(Interpolate, data, fm) {
  this->utils.setFieldData(bf, fm);
  this->utils.setFieldData(gbf, fm);
  for (int i = 0; i < num_fields; ++i) {
    this->utils.setFieldData(nodal[i], fm);
    this->utils.setFieldData(u[i], fm);
    this->utils.setFieldData(gu[i], fm);
  }
  (void)data;
}

PHX_EVALUATE_FIELDS(Interpolate, workset) {
  for (int elem = 0; elem < workset.size; ++elem) {
    for (int ip = 0; ip < num_ips; ++ip) {
      for (int f = 0; f < num_fields; ++f) {
        u(elem, ip) = nodal(elem, 0) * bf(elem, 0, ip);
        for (int n = 1; n < num_nodes; ++n)
          u(elem, ip) += nodal(elem, n) * bf(elem, n, ip);
        for (int i = 0; i < num_dims; ++i) {
          gu(elem, ip, i) = nodal(elem, 0) * gbf(elem, 0, ip, i);
          for (int n = 1; n < num_nodes; ++n)
            gu(elem, n, i) += nodal(elem, n) * gbf(elem, n, ip, i);
  }}}}
}

} // end namespace goal
