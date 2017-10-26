#include "goal_disc.hpp"
#include "goal_sol_info.hpp"

namespace goal {

SolInfo::SolInfo(Disc* d) {
  disc = d;
  owned = new LinearObj;
  ghost = new LinearObj;
  auto owned_map = disc->get_owned_map();
  auto ghost_map = disc->get_ghost_map();
  auto owned_graph = disc->get_owned_graph();
  auto ghost_graph = disc->get_ghost_graph();
  importer = rcp(new ImportT(owned_map, ghost_map));
  exporter = rcp(new ExportT(ghost_map, owned_map));
  owned->R = rcp(new VectorT(owned_map));
  owned->dMdu = rcp(new VectorT(owned_map));
  owned->dRdu = rcp(new MatrixT(owned_graph));
  ghost->R = rcp(new VectorT(ghost_map));
  ghost->dMdu = rcp(new VectorT(ghost_map));
  ghost->dRdu = rcp(new MatrixT(ghost_graph));
}

SolInfo::~SolInfo() {
  delete ghost;
  delete owned;
}

Disc* SolInfo::get_disc() {
  return disc;
}

void SolInfo::gather_R() {
  owned->R->doExport(*(ghost->R), *exporter, Tpetra::ADD);
}

void SolInfo::gather_dMdu() {
  owned->dMdu->doExport(*(ghost->dMdu), *exporter, Tpetra::ADD);
}

void SolInfo::gather_dRdu() {
  owned->dRdu->doExport(*(ghost->dRdu), *exporter, Tpetra::ADD);
}

void SolInfo::gather_all() {
  gather_R();
  gather_dMdu();
  gather_dRdu();
}

void SolInfo::zero_R() {
  owned->R->putScalar(0.0);
  ghost->R->putScalar(0.0);
}

void SolInfo::zero_dMdu() {
  owned->dMdu->putScalar(0.0);
  ghost->dMdu->putScalar(0.0);
}

void SolInfo::zero_dRdu() {
  owned->dRdu->setAllToScalar(0.0);
  ghost->dRdu->setAllToScalar(0.0);
}

void SolInfo::zero_all() {
  zero_R();
  zero_dRdu();
  zero_dMdu();
}

void SolInfo::resume_fill() {
  owned->dRdu->resumeFill();
  ghost->dRdu->resumeFill();
}

void SolInfo::complete_fill() {
  owned->dRdu->fillComplete();
  ghost->dRdu->fillComplete();
}

SolInfo* create_sol_info(Disc* d) {
  return new SolInfo(d);
}

void destroy_sol_info(SolInfo* s) {
  delete s;
}

}
