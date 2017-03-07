#include <apfNumbering.h>

#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_solution_info.hpp"

namespace goal {

using Teuchos::rcp;

SolutionInfo::SolutionInfo(RCP<Indexer> indexer, int n) {
  num_qois = n;
  owned = rcp(new LinearObj);
  ghost = rcp(new LinearObj);
  auto om = indexer->get_owned_map();
  auto gm = indexer->get_ghost_map();
  auto og = indexer->get_owned_graph();
  auto gg = indexer->get_ghost_graph();
  importer = rcp(new Import(om, gm));
  exporter = rcp(new Export(gm, om));
  owned->u = rcp(new Vector(om));
  owned->R = rcp(new Vector(om));
  owned->dJdu = rcp(new Vector(om, num_qois));
  owned->z = rcp(new Vector(om, num_qois));
  owned->dRdu = rcp(new Matrix(og));
  ghost->u = rcp(new Vector(gm));
  ghost->R = rcp(new Vector(gm));
  ghost->dRdu = rcp(new Matrix(gg));
  ghost->dJdu = rcp(new Vector(gm, num_qois));
  ghost->z = rcp(new Vector(gm, num_qois));
}

void SolutionInfo::gather_u() {
  owned->u->doExport(*(ghost->u), *exporter, Tpetra::INSERT);
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

void SolutionInfo::scatter_u() {
  ghost->u->doImport(*(owned->u), *importer, Tpetra::INSERT);
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

static void fill_field(RCP<Field> f, RCP<Indexer> i, RCP<Vector> x) {
  auto idx = f->get_associated_dof_idx();
  auto apf_field = f->get_apf_field();
  auto components = f->get_num_components();
  auto numbering = i->get_apf_numbering(idx);
  auto data = x->get1dView();
  apf::DynamicArray<apf::Node> owned;
  apf::getNodes(numbering, owned);
  double values[3] = {0.0, 0.0, 0.0};
  for (size_t n = 0; n < owned.size(); ++n) {
    auto node = owned[n];
    auto ent = node.entity;
    auto entn = node.node;
    for (int c = 0; c < components; ++c) {
      LO row = i->get_owned_lid(idx, node, c);
      values[c] = data[row];
    }
    apf::setComponents(apf_field, ent, entn, values);
  }
  apf::synchronize(apf_field);
}

void fill_fields(
    std::vector<RCP<Field> > fields, RCP<Indexer> indexer, RCP<Vector> x) {
  double t0 = time();
  for (size_t i = 0; i < fields.size(); ++i) fill_field(fields[i], indexer, x);
  double t1 = time();
  print(" > fields filled in %f seconds", t1 - t0);
}

}  // namespace goal
