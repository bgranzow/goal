#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_physics.hpp"

namespace goal {

using Teuchos::rcp;

Physics::Physics(Discretization* d) {
  disc = d;
}

Physics::~Physics() {
}

void Physics::build_coarse_indexer() {
  indexer = create_indexer(disc, u);
}

void Physics::build_fine_indexer() {
  indexer = create_indexer(disc, u_fine);
}

void Physics::build_error_indexer() {
  indexer = create_indexer(disc, e);
}

void Physics::destroy_indexer() {
  goal::destroy_indexer(indexer);
}

void Physics::build_primal_model() {
  auto t0 = time();
  auto nes = disc->get_num_elem_sets();
  auto nss = disc->get_num_side_sets();
  vfms.resize(nes);
  nfms.resize(nss);
  for (int i = 0; i < nes; ++i) {
    elem_set = i;
    vfms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_primal_volumetric(vfms[i]);
  }
  for (int i = 0; i < nss; ++i) {
    side_set = i;
    nfms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_primal_neumann(nfms[i]);
  }
  dfm = rcp(new PHX::FieldManager<goal::Traits>);
  build_primal_dirichlet(dfm);
  auto t1 = time();
  print(" > primal model built in %f seconds", t1 - t0);
}

void Physics::build_dual_model() {
  auto t0 = time();
  auto nes = disc->get_num_elem_sets();
  auto nss = disc->get_num_side_sets();
  vfms.resize(nes);
  nfms.resize(nss);
  for (int i = 0; i < nes; ++i) {
    elem_set = i;
    vfms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_dual_volumetric(vfms[i]);
  }
  for (int i = 0; i < nss; ++i) {
    side_set = i;
    nfms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_dual_neumann(nfms[i]);
  }
  dfm = rcp(new PHX::FieldManager<goal::Traits>);
  build_dual_dirichlet(dfm);
  auto t1 = time();
  print(" > dual model built in %f seconds", t1 - t0);
}

void Physics::build_error_model() {
  auto t0 = time();
  auto nes = disc->get_num_elem_sets();
  auto nss = disc->get_num_side_sets();
  vfms.resize(nes);
  nfms.resize(nss);
  for (int i = 0; i < nes; ++i) {
    elem_set = i;
    vfms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_error_volumetric(vfms[i]);
  }
  for (int i = 0; i < nss; ++i) {
    side_set = i;
    nfms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_error_neumann(nfms[i]);
  }
  dfm = rcp(new PHX::FieldManager<goal::Traits>);
  build_error_dirichlet(dfm);
  auto t1 = time();
  print(" > error model built in %f seconds", t1 - t0);
}

void Physics::destroy_model() {
  for (int i = 0; i < vfms.size(); ++i)
    vfms[i] = Teuchos::null;
  for (int i = 0; i < nfms.size(); ++i)
    nfms[i] = Teuchos::null;
  dfm = Teuchos::null;
}

static int get_q_degree(const int p) {
  int q = 0;
  switch (p) {
    case 1: q = 1; break;
    case 2: q = 2; break;
    case 3: q = 4; break;
  }
  return q;
}

static int get_highest_p(std::vector<goal::Field*> const& u) {
  int p = 0;
  for (size_t i = 0; i < u.size(); ++i)
    p = std::max(p, u[i]->get_p_order());
  return p;
}

void Physics::build_enriched_primal_fields() {
  GOAL_DEBUG_ASSERT(u_fine.size() == 0);
  int highest_p = get_highest_p(u);
  int q = get_q_degree(highest_p + 1);
  for (size_t i = 0; i < u.size(); ++i) {
    auto idx = u[i]->get_associated_dof_idx();
    auto n = u[i]->name() + "_fine";
    auto b = u[i]->get_basis_type();
    auto p = u[i]->get_p_order() + 1;
    FieldInfo info = {disc, n, p, q, b};
    u_fine.push_back(goal::create_field(info));
    u_fine.back()->set_associated_dof_idx(idx);
    u_fine.back()->set_seed(1.0);
    project_field(u_fine.back(), u[i]);
  }
}

void Physics::build_dual_fields() {
  GOAL_DEBUG_ASSERT(z.size() == 0);
  for (size_t i = 0; i < u.size(); ++i) {
    auto idx = u[i]->get_associated_dof_idx();
    auto n = u[i]->name() + "_dual";
    auto b = u[i]->get_basis_type();
    auto p = u[i]->get_p_order();
    auto q = u[i]->get_q_degree();
    FieldInfo info = {disc, n, p, q, b};
    z.push_back(goal::create_field(info));
    z.back()->set_associated_dof_idx(idx);
    z.back()->set_seed(0.0);
  }
}

void Physics::build_enriched_dual_fields() {
  GOAL_DEBUG_ASSERT(z_fine.size() == 0);
  int highest_p = get_highest_p(u);
  int q = get_q_degree(highest_p + 1);
  for (size_t i = 0; i < u.size(); ++i) {
    auto idx = u[i]->get_associated_dof_idx();
    auto n = u[i]->name() + "_dual_fine";
    auto b = u[i]->get_basis_type();
    auto p = u[i]->get_p_order() + 1;
    FieldInfo info = {disc, n, p, q, b};
    z_fine.push_back(goal::create_field(info));
    z_fine.back()->set_associated_dof_idx(idx);
    z_fine.back()->set_seed(0.0);
    project_field(z_fine.back(), z[i]);
  }
}

void Physics::build_error_fields() {
  GOAL_DEBUG_ASSERT(e.size() == 0);
  for (size_t i = 0; i < u.size(); ++i) {
    auto idx = u[i]->get_associated_dof_idx();
    auto n = u[i]->name() + "_error";
    FieldInfo info = {disc, n, 1, 1, LAGRANGE};
    e.push_back(goal::create_field(info));
    e.back()->set_associated_dof_idx(idx);
    e.back()->set_seed(0.0);
  }
}

void Physics::destroy_enriched_primal_fields() {
  for (size_t i = 0; i < u_fine.size(); ++i)
    destroy_field(u_fine[i]);
  u_fine.resize(0);
}

void Physics::destroy_dual_fields() {
  for (size_t i = 0; i < z.size(); ++i)
    destroy_field(z[i]);
  z.resize(0);
}

void Physics::destroy_enriched_dual_fields() {
  for (size_t i = 0; i < z_fine.size(); ++i)
    destroy_field(z_fine[i]);
  z_fine.resize(0);
}

void Physics::destroy_error_fields() {
  for (size_t i = 0; i < e.size(); ++i)
    destroy_field(e[i]);
  e.resize(0);
}

void set_extended_data_type_dims(Indexer* indexer, FieldManager fm, int t) {
  std::vector<PHX::index_size_type> dd;
  dd.push_back(indexer->get_num_total_dofs(t));
  fm->setKokkosExtendedDataTypeDimensions<goal::Traits::Jacobian>(dd);
}

} // end namespace goal
