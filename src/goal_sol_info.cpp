#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_sol_info.hpp"

namespace goal {

using Teuchos::rcp;

SolInfo::SolInfo(Indexer* i, int nqoi) {
  owned = new LinearObj;
  ghost = new LinearObj;
  auto om = i->get_owned_map();
  auto gm = i->get_ghost_map(); 
  auto og = i->get_owned_graph();
  auto gg = i->get_ghost_graph();
  importer = rcp(new Import(om, gm));
  exporter = rcp(new Export(gm, om));
  owned->du = rcp(new Vector(om));
  ghost->du = rcp(new Vector(gm));
  owned->R = rcp(new Vector(om));
  ghost->R = rcp(new Vector(gm));
  owned->dRdu = rcp(new Matrix(og));
  ghost->dRdu = rcp(new Matrix(gg));
  if (nqoi > 0) {
    owned->dJdu = rcp(new MultiVector(om, nqoi));
    ghost->dJdu = rcp(new MultiVector(gm, nqoi));
  }
}

SolInfo::~SolInfo() {
  delete owned;
  delete ghost;
}

void SolInfo::gather_du() {
  owned->du->doExport(*(ghost->du), *exporter, Tpetra::INSERT);
}

void SolInfo::gather_R() {
  owned->R->doExport(*(ghost->R), *exporter, Tpetra::ADD);
}

void SolInfo::gather_dRdu() {
  owned->dRdu->doExport(*(ghost->dRdu), *exporter, Tpetra::ADD);
}

void SolInfo::gather_dJdu() {
  GOAL_DEBUG_ASSERT( Teuchos::nonnull(owned->dJdu) );
  GOAL_DEBUG_ASSERT( Teuchos::nonnull(ghost->dJdu) );
  owned->dJdu->doExport(*(ghost->dJdu), *exporter, Tpetra::ADD);
}

void SolInfo::scatter_du() {
  ghost->du->doImport(*(owned->du), *importer, Tpetra::INSERT);
}

void SolInfo::scatter_R() {
  ghost->R->doImport(*(owned->R), *importer, Tpetra::INSERT);
}

void SolInfo::scatter_dRdu() {
  ghost->dRdu->doImport(*(owned->R), *importer, Tpetra::INSERT);
}

void SolInfo::scatter_dJdu() {
  GOAL_DEBUG_ASSERT( Teuchos::nonnull(owned->dJdu) );
  GOAL_DEBUG_ASSERT( Teuchos::nonnull(ghost->dJdu) );
  ghost->dJdu->doImport(*(owned->dJdu), *importer, Tpetra::INSERT);
}

SolInfo* create_sol_info(Indexer* i, int nqoi) {
  return new SolInfo(i, nqoi);
}

void destroy_sol_info(SolInfo* s) {
  delete s;
}

} // end namespace goal
