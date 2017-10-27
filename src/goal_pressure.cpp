#include <apfMesh.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_eval_modes.hpp"
#include "goal_pressure.hpp"
#include "goal_sol_info.hpp"

namespace goal {

Pressure<ST>::Pressure(apf::Field* base, int mode) {
  disc = 0;
  elem = 0;
  num_dims = 0;
  num_nodes = 0;
  field = base;
  shape = apf::getShape(field);
  num_dims = apf::getMesh(field)->getDimension();
  if (mode == PRIMAL) op = &Pressure<ST>::scatter_primal;
  else if (mode == NONE) op = &Pressure<ST>::scatter_none;
  else fail("displacement: invalid mode: %d", mode);
  auto fname = (std::string)apf::getName(base);
  this->name = fname.substr(0, 1);
}

Pressure<ST>::~Pressure() {
}

ST& Pressure<ST>::val() {
  return value;
}

ST& Pressure<ST>::grad(int i) {
  return gradient[i];
}

ST& Pressure<ST>::nodal(int n) {
  return node[n];
}

ST& Pressure<ST>::resid(int n) {
  return residual[n];
}

void Pressure<ST>::pre_process(SolInfo* s) {
  disc = s->get_disc();
}

void Pressure<ST>::gather(apf::MeshElement* me) {
  auto ent = apf::getMeshEntity(me);
  num_nodes = disc->get_num_nodes(ent);
  residual.resize(num_nodes);
  elem = apf::createElement(field, me);
  apf::getScalarNodes(elem, node);
  for (int n = 0; n < num_nodes; ++n)
    residual[n] = 0.0;
}

void Pressure<ST>::at_point(apf::Vector3 const& p, double, double) {
  auto me = apf::getMeshElement(elem);
  apf::getBF(shape, me, p, BF);
  apf::getGradBF(shape, me, p, GBF);
  value = nodal(0) * BF[0];
  for (int n = 1; n < num_nodes; ++n)
    value += nodal(n) * BF[n];
  for (int i = 0; i < num_dims; ++i) {
    gradient[i] = nodal(0) * GBF[0][i];
    for (int n = 1; n < num_nodes; ++n)
      gradient[i] += nodal(n) * GBF[n][i];
  }
}

void Pressure<ST>::scatter_none(SolInfo*) {
}

void Pressure<ST>::scatter_primal(SolInfo* s) {
  auto ent = apf::getMeshEntity(elem);
  auto R = s->ghost->R->get1dViewNonConst();
  for (int n = 0; n < num_nodes; ++n) {
    LO row = disc->get_lid(ent, n, num_dims);
    R[row] += resid(n);
  }
}

void Pressure<ST>::scatter(SolInfo* s) {
  op(this, s);
  apf::destroyElement(elem);
  elem = 0;
}

void Pressure<ST>::post_process(SolInfo*) {
  disc = 0;
}

Pressure<FADT>::Pressure(apf::Field* base, int mode) {
  disc = 0;
  elem = 0;
  num_dims = 0;
  num_nodes = 0;
  num_dofs = 0;
  field = base;
  shape = apf::getShape(field);
  num_dims = apf::getMesh(field)->getDimension();
  if (mode == NONE) op = &Pressure<FADT>::scatter_none;
  else if (mode == PRIMAL) op = &Pressure<FADT>::scatter_primal;
  else if (mode == ADJOINT) op = &Pressure<FADT>::scatter_adjoint;
  else fail("displacement: invalid mode: %d", mode);
  auto fname = (std::string)apf::getName(base);
  this->name = fname.substr(0, 1);
}

Pressure<FADT>::~Pressure() {
}

FADT& Pressure<FADT>::val() {
  return value;
}

FADT& Pressure<FADT>::grad(int i) {
  return gradient[i];
}

FADT& Pressure<FADT>::nodal(int n) {
  return node_fadt[n];
}

FADT& Pressure<FADT>::resid(int n) {
  return residual[n];
}

void Pressure<FADT>::pre_process(SolInfo* s) {
  disc = s->get_disc();
  gradient.resize(num_dims);
}

void Pressure<FADT>::gather(apf::MeshElement* me) {
  auto ent = apf::getMeshEntity(me);
  int num_eqs = disc->get_num_eqs();
  num_nodes = disc->get_num_nodes(ent);
  num_dofs = disc->get_num_dofs(ent);
  node_fadt.resize(num_nodes);
  residual.resize(num_nodes);
  elem = apf::createElement(field, me);
  apf::getScalarNodes(elem, node_st);
  for (int n = 0; n < num_nodes; ++n) {
    int eq = n*num_eqs + num_dims;
    nodal(n).diff(eq, num_dofs);
    nodal(n).val() = node_st[n];
    resid(n) = 0.0;
  }
}

void Pressure<FADT>::at_point(apf::Vector3 const& p, double, double) {
  auto me = apf::getMeshElement(elem);
  apf::getBF(shape, me, p, BF);
  apf::getGradBF(shape, me, p, GBF);
  value = nodal(0) * BF[0];
  for (int n = 1; n < num_nodes; ++n)
    value += nodal(n) * BF[n];
  for (int i = 0; i < num_dims; ++i) {
    gradient[i] = nodal(0) * GBF[0][i];
    for (int n = 1; n < num_nodes; ++n)
      gradient[i] += nodal(n) * GBF[n][i];
  }
}

void Pressure<FADT>::scatter_none(SolInfo*) {
}

void Pressure<FADT>::scatter_primal(SolInfo* s) {
  using Teuchos::arrayView;
  auto ent = apf::getMeshEntity(elem);
  auto R = s->ghost->R->get1dViewNonConst();
  auto dRdu = s->ghost->dRdu;
  std::vector<LO> cols(num_dofs);
  disc->get_lids(ent, cols);
  auto c = arrayView(&cols[0], num_dofs);
  for (int n = 0; n < num_nodes; ++n) {
    auto v = resid(n);
    auto view = arrayView(&(v.fastAccessDx(0)), num_dofs);
    LO row = disc->get_lid(ent, n, num_dims);
    R[row] += v.val();
    dRdu->sumIntoLocalValues(row, c, view, num_dofs);
  }
}

void Pressure<FADT>::scatter_adjoint(SolInfo* s) {
  using Teuchos::arrayView;
  auto ent = apf::getMeshEntity(elem);
  auto R = s->ghost->R->get1dViewNonConst();
  auto dRduT = s->ghost->dRdu;
  std::vector<LO> cols(num_dofs);
  disc->get_lids(ent, cols);
  for (int n = 0; n < num_nodes; ++n) {
    auto v = resid(n);
    auto view = arrayView(&(v.fastAccessDx(0)), num_dofs);
    LO row = disc->get_lid(ent, n, num_dims);
    R[row] += v.val();
    for (int dof = 0; dof < num_dofs; ++dof)
      dRduT->sumIntoLocalValues(
          cols[dof], arrayView(&row, 1), arrayView(&view[dof], 1));
  }
}

void Pressure<FADT>::scatter(SolInfo* s) {
  op(this, s);
  apf::destroyElement(elem);
  elem = 0;
}

void Pressure<FADT>::post_process(SolInfo*) {
  disc = 0;
}

}
