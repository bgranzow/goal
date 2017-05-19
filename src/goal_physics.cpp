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

void set_extended_data_type_dims(Indexer* indexer, FieldManager fm, int t) {
  std::vector<PHX::index_size_type> dd;
  dd.push_back(indexer->get_num_total_dofs(t));
  fm->setKokkosExtendedDataTypeDimensions<goal::Traits::Jacobian>(dd);
}

} // end namespace goal
