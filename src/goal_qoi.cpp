#include <apf.h>
#include <PCU.h>

#include "goal_disc.hpp"
#include "goal_qoi.hpp"
#include "goal_sol_info.hpp"

namespace goal {

QoI<ST>::QoI() :
    elem(0),
    disc(0),
    qoi_value(0.0),
    elem_value(0.0) {
}

QoI<ST>::~QoI() {
}

void QoI<ST>::pre_process(SolInfo* s) {
  qoi_value = 0.0;
  disc = s->get_disc();
}

void QoI<ST>::gather(apf::MeshElement* me) {
  elem = me;
  elem_value = 0.0;
}

void QoI<ST>::scatter(SolInfo*) {
  qoi_value += elem_value;
  elem = 0;
}

void QoI<ST>::post_process(SolInfo*) {
  disc = 0;
  PCU_Add_Doubles(&qoi_value, 1);
}

QoI<FADT>::QoI() :
    elem(0),
    disc(0),
    qoi_value(0.0),
    elem_value(0.0) {
}

QoI<FADT>::~QoI() {
}

void QoI<FADT>::pre_process(SolInfo* s) {
  qoi_value = 0.0;
  disc = s->get_disc();
}

void QoI<FADT>::gather(apf::MeshElement* me) {
  elem = me;
  auto ent = apf::getMeshEntity(elem);
  auto num_dofs = disc->get_num_dofs(ent);
  elem_value.diff(0, num_dofs);
  elem_value.fastAccessDx(0) = 0.0;
}

void QoI<FADT>::scatter(SolInfo* s) {
  std::vector<LO> rows;
  auto dMdu = s->ghost->dMdu->get1dViewNonConst();
  auto ent = apf::getMeshEntity(elem);
  auto num_dofs = disc->get_num_dofs(ent);
  disc->get_lids(ent, rows);
  for (int dof = 0; dof < num_dofs; ++dof) {
    LO row = rows[dof];
    auto val = get_elem_value().fastAccessDx(dof);
    dMdu[row] += val;
  }
  qoi_value += elem_value.val();
  elem = 0;
}

void QoI<FADT>::post_process(SolInfo*) {
  disc = 0;
  PCU_Add_Doubles(&qoi_value, 1);
}

}
