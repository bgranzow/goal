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

static void add_to_field(RCP<Field> f, RCP<Indexer> i, RCP<Vector> du) {
  auto idx = f->get_associated_dof_idx();
  auto apf_field = f->get_apf_field();
  auto components = f->get_num_components();
  auto numbering = i->get_apf_numbering(idx);
  auto data = du->get1dView();
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(numbering, owned);
  double values[3] = {0.0, 0.0, 0.0};
  for (size_t n = 0; n < owned.size(); ++n) {
    auto node = owned[n];
    auto ent = node.entity;
    auto entn = node.node;
    apf::getComponents(apf_field, ent, entn, values);
    for (int c = 0; c < components; ++c) {
      LO row = i->get_owned_lid(idx, node, c);
      values[c] += data[row];
    }
    apf::setComponents(apf_field, ent, entn, values);
  }
  apf::synchronize(apf_field);
}

void add_to_fields(
    std::vector<RCP<Field> > const& f, RCP<Indexer> i, RCP<Vector> du) {
  double t0 = time();
  for (size_t j = 0; j < f.size(); ++j)
    add_to_field(f[j], i, du);
  double t1 = time();
  print(" > added to fields in %f seconds", t1 - t0);
}

}  /* namespace goal */
