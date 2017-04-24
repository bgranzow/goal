#include <apfNumbering.h>

#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_solution_info.hpp"

namespace goal {

using Teuchos::rcp;

SolutionInfo::SolutionInfo(RCP<Indexer> indexer) {
  owned = rcp(new LinearObj);
  ghost = rcp(new LinearObj);
  auto om = indexer->get_owned_map();
  auto gm = indexer->get_ghost_map();
  auto og = indexer->get_owned_graph();
  auto gg = indexer->get_ghost_graph();
  importer = rcp(new Import(om, gm));
  exporter = rcp(new Export(gm, om));
  owned->du = rcp(new Vector(om));
  owned->R = rcp(new Vector(om));
  owned->dJdu = rcp(new Vector(om));
  owned->z = rcp(new Vector(om));
  owned->dRdu = rcp(new Matrix(og));
  ghost->du = rcp(new Vector(gm));
  ghost->R = rcp(new Vector(gm));
  ghost->dRdu = rcp(new Matrix(gg));
  ghost->dJdu = rcp(new Vector(gm));
  ghost->z = rcp(new Vector(gm));
}

void SolutionInfo::gather_du() {
  owned->du->doExport(*(ghost->du), *exporter, Tpetra::INSERT);
}

void SolutionInfo::gather_R() {
  owned->R->doExport(*(ghost->R), *exporter, Tpetra::ADD);
}

void SolutionInfo::gather_dRdu() {
  owned->dRdu->doExport(*(ghost->dRdu), *exporter, Tpetra::ADD);
}

void SolutionInfo::gather_dJdu() {
  owned->dJdu->doExport(*(ghost->dJdu), *exporter, Tpetra::ADD);
}

void SolutionInfo::gather_z() {
  owned->z->doExport(*(ghost->z), *exporter, Tpetra::ADD);
}

void SolutionInfo::scatter_du() {
  ghost->du->doImport(*(owned->du), *importer, Tpetra::INSERT);
}

void SolutionInfo::scatter_R() {
  ghost->R->doImport(*(owned->R), *importer, Tpetra::INSERT);
}

void SolutionInfo::scatter_dRdu() {
  ghost->dRdu->doImport(*(owned->R), *importer, Tpetra::INSERT);
}

void SolutionInfo::scatter_dJdu() {
  ghost->dJdu->doImport(*(owned->dJdu), *importer, Tpetra::INSERT);
}

void SolutionInfo::scatter_z() {
  ghost->z->doImport(*(owned->z), *importer, Tpetra::INSERT);
}

}  /* namespace goal */
