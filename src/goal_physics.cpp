#include "goal_physics.hpp"
#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"

namespace goal {

using Teuchos::rcp;

Physics::Physics(RCP<Discretization> d) { disc = d; }

Physics::~Physics() {}

void Physics::build_coarse_indexer() { indexer = rcp(new Indexer(disc, u)); }

void Physics::build_fine_indexer() { indexer = rcp(new Indexer(disc, u_fine)); }

void Physics::destroy_indexer() { indexer = Teuchos::null; }

void Physics::build_primal_model() {
  auto t0 = time();
  auto neb = disc->get_num_elem_blocks();
  volumetric_fms.resize(neb);
  for (int i = 0; i < neb; ++i) {
    indexer->set_elem_block(i);
    volumetric_fms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_primal_volumetric(volumetric_fms[i]);
  }
  neumann_fm = rcp(new PHX::FieldManager<goal::Traits>);
  build_primal_neumann(neumann_fm);
  dirichlet_fm = rcp(new PHX::FieldManager<goal::Traits>);
  build_primal_dirichlet(dirichlet_fm);
  auto t1 = time();
  print(" > primal model built in %f seconds", t1 - t0);
}

void Physics::build_dual_model() {
  auto t0 = time();
  auto neb = disc->get_num_elem_blocks();
  volumetric_fms.resize(neb);
  for (int i = 0; i < neb; ++i) {
    indexer->set_elem_block(i);
    volumetric_fms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_dual_volumetric(volumetric_fms[i]);
  }
  neumann_fm = rcp(new PHX::FieldManager<goal::Traits>);
  build_dual_neumann(neumann_fm);
  dirichlet_fm = rcp(new PHX::FieldManager<goal::Traits>);
  build_dual_dirichlet(dirichlet_fm);
  auto t1 = time();
  print(" > dual model built in %f seconds", t1 - t0);
}

void Physics::build_error_model() {
  auto t0 = time();
  auto neb = disc->get_num_elem_blocks();
  volumetric_fms.resize(neb);
  for (int i = 0; i < neb; ++i) {
    indexer->set_elem_block(i);
    volumetric_fms[i] = rcp(new PHX::FieldManager<goal::Traits>);
    build_error_volumetric(volumetric_fms[i]);
  }
  neumann_fm = rcp(new PHX::FieldManager<goal::Traits>);
  build_error_neumann(neumann_fm);
  dirichlet_fm = rcp(new PHX::FieldManager<goal::Traits>);
  build_error_dirichlet(dirichlet_fm);
  auto t1 = time();
  print(" > error model built in %f seconds", t1 - t0);
}

void Physics::destroy_model() {
  for (int i = 0; i < volumetric_fms.size(); ++i)
    volumetric_fms[i] = Teuchos::null;
  neumann_fm = Teuchos::null;
  dirichlet_fm = Teuchos::null;
}

static void build_enriched_primal_field(int q, RCP<Field> u,
    RCP<Discretization> disc, std::vector<RCP<Field> >& u_fine) {
  auto idx = u->get_associated_dof_idx();
  auto n = u->get_name() + "_fine";
  auto t = u->get_value_type();
  auto b = u->get_basis_type();
  auto p = u->get_p_order() + 1;
  FieldInfo i = {disc, n, p, q, t, b};
  u_fine.push_back(rcp(new Field(&i)));
  u_fine.back()->set_associated_dof_idx(idx);
  u_fine.back()->set_dof_status(true);
  u_fine.back()->set_seed(1.0);
  project_field(u_fine.back(), u);
}

static void build_dual_field(int q, RCP<Field> u, RCP<Discretization> disc,
    std::vector<RCP<Field> >& z) {
  auto idx = u->get_associated_dof_idx();
  auto n = u->get_name() + "_dual";
  auto t = u->get_value_type();
  auto b = u->get_basis_type();
  auto p = u->get_p_order();
  FieldInfo i = {disc, n, p, q, t, b};
  z.push_back(rcp(new Field(&i)));
  z.back()->set_associated_dof_idx(idx);
  z.back()->set_dof_status(false);
  z.back()->set_seed(0.0);
}

static void build_enriched_dual_field(int q, RCP<Field> u,
    RCP<Discretization> disc, std::vector<RCP<Field> >& z_fine) {
  auto idx = u->get_associated_dof_idx();
  auto n = u->get_name() + "_dual_fine";
  auto t = u->get_value_type();
  auto b = u->get_basis_type();
  auto p = u->get_p_order() + 1;
  FieldInfo i = {disc, n, p, q, t, b};
  z_fine.push_back(rcp(new Field(&i)));
  z_fine.back()->set_associated_dof_idx(idx);
  z_fine.back()->set_dof_status(false);
  z_fine.back()->set_seed(0.0);
}

static void build_error_field(
    RCP<Field> u, RCP<Discretization> disc, std::vector<RCP<Field> >& e) {
  auto idx = u->get_associated_dof_idx();
  auto n = u->get_name() + "_error";
  auto t = u->get_value_type();
  FieldInfo i = {disc, n, 1, 1, t, LAGRANGE};
  e.push_back(rcp(new Field(&i)));
  e.back()->set_associated_dof_idx(idx);
  e.back()->set_dof_status(false);
  e.back()->set_seed(0.0);
}

static int get_q_degree(int p) {
  int q = 0;
  switch (p) {
    case 1:
      q = 1;
      break;
    case 2:
      q = 2;
      break;
    case 3:
      q = 4;
      break;
  }
  return q;
}

void Physics::build_enriched_data() {
  assert(u_fine.size() == 0);
  assert(z.size() == 0);
  assert(z_fine.size() == 0);
  assert(e.size() == 0);
  int highest_p = 0;
  for (size_t i = 0; i < u.size(); ++i)
    highest_p = std::max(highest_p, u[i]->get_p_order());
  int q = get_q_degree(highest_p + 1);
  for (size_t i = 0; i < u.size(); ++i) {
    build_enriched_primal_field(q, u[i], disc, u_fine);
    build_dual_field(q, u[i], disc, z);
    build_enriched_dual_field(q, u[i], disc, z_fine);
    build_error_field(u[i], disc, e);
  }
}

void Physics::destroy_enriched_data() {
  for (size_t i = 0; i < u_fine.size(); ++i) u_fine[i] = Teuchos::null;
  for (size_t i = 0; i < z.size(); ++i) z[i] = Teuchos::null;
  for (size_t i = 0; i < z_fine.size(); ++i) z_fine[i] = Teuchos::null;
  for (size_t i = 0; i < e.size(); ++i) e[i] = Teuchos::null;
  u_fine.resize(0);
  z.resize(0);
  z_fine.resize(0);
  e.resize(0);
}

void set_extended_data_type_dims(RCP<Indexer> indexer, FieldManager fm) {
  {
    std::vector<PHX::index_size_type> dd;
    dd.push_back(indexer->get_num_total_elem_dofs());
    fm->setKokkosExtendedDataTypeDimensions<goal::Traits::Jacobian>(dd);
  }
}

}  // namespace goal
