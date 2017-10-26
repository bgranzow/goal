#include <apfMesh.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_eval_modes.hpp"
#include "goal_displacement.hpp"
#include "goal_sol_info.hpp"

namespace goal {

Displacement<ST>::Displacement(apf::Field* base, const int mode) {
  disc = 0;
  elem = 0;
  num_nodes = 0;
  field = base;
  shape = apf::getShape(field);
  num_dims = apf::getMesh(field)->getDimension();
  if (mode == PRIMAL) op = &Displacement<ST>::scatter_primal;
  else if (mode == NONE) op = &Displacement<ST>::scatter_none;
  else fail("displacement: invalid mode: %d", mode);
  auto fname = (std::string)apf::getName(base);
  this->name = fname.substr(0, 1);
}

Displacement<ST>::~Displacement() {
}

ST& Displacement<ST>::val(const int i) {
  return value[i];
}

ST& Displacement<ST>::grad(const int i, const int j) {
  return gradient[i][j];
}

ST& Displacement<ST>::nodal(const int n, const int i) {
  return node[n][i];
}

ST& Displacement<ST>::resid(const int n, const int i) {
  return residual[n * num_dims + i];
}

void Displacement<ST>::pre_process(SolInfo* s) {
  disc = s->get_disc();
}

void Displacement<ST>::gather(apf::MeshElement* me) {
  auto ent = apf::getMeshEntity(me);
  num_nodes = disc->get_num_nodes(ent);
  elem = apf::createElement(field, me);
  apf::getVectorNodes(elem, node);
  residual.resize(num_dims * num_nodes);
  for (int i = 0; i < num_dims * num_nodes; ++i)
    residual[i] = 0.0;
}

void Displacement<ST>::at_point(apf::Vector3 const& p, double, double) {
  auto me = apf::getMeshElement(elem);
  apf::getBF(shape, me, p, BF);
  apf::getGradBF(shape, me, p, GBF);
  for (int i = 0; i < num_dims; ++i) {
    value[i] = nodal(0, i) * BF[0];
    for (int n = 1; n < num_nodes; ++n)
      value[i] += nodal(n, i) * BF[n];
    for (int j = 0; j < num_dims; ++j) {
      gradient[i][j] = nodal(0, i) * GBF[0][j];
      for (int n = 1; n < num_nodes; ++n)
        gradient[i][j] += nodal(n, i) * GBF[n][j];
    }
  }
}

void Displacement<ST>::scatter_none(SolInfo*) {
}

void Displacement<ST>::scatter_primal(SolInfo* s) {
  auto ent = apf::getMeshEntity(elem);
  auto R = s->ghost->R->get1dViewNonConst();
  for (int n = 0; n < num_nodes; ++n) {
    for (int d = 0; d < num_dims; ++d) {
      LO row = disc->get_lid(ent, n, d);
      R[row] += resid(n, d);
    }
  }
}

void Displacement<ST>::scatter(SolInfo* s) {
  op(this, s);
  apf::destroyElement(elem);
  elem = 0;
}

void Displacement<ST>::post_process(SolInfo*) {
  disc = 0;
}

Displacement<FADT>::Displacement(apf::Field* base, const int mode) {
  disc = 0;
  elem = 0;
  num_nodes = 0;
  num_dofs = 0;
  field = base;
  shape = apf::getShape(field);
  num_dims = apf::getMesh(field)->getDimension();
  if (mode == NONE) op = &Displacement<FADT>::scatter_none;
  else if (mode == PRIMAL) op = &Displacement<FADT>::scatter_primal;
  else if (mode == ADJOINT) op = &Displacement<FADT>::scatter_adjoint;
  else fail("displacement: invalid mode: %d", mode);
  auto fname = (std::string)apf::getName(base);
  this->name = fname.substr(0, 1);
}

Displacement<FADT>::~Displacement() {
}

FADT& Displacement<FADT>::val(const int i) {
  return value[i];
}

FADT& Displacement<FADT>::grad(const int i, const int j) {
  return gradient[i * num_dims + j];
}

FADT& Displacement<FADT>::nodal(const int n, const int i) {
  return node_fadt[n * num_dims + i];
}

FADT& Displacement<FADT>::resid(const int n, const int i) {
  return residual[n * num_dims + i];
}

void Displacement<FADT>::pre_process(SolInfo* s) {
  disc = s->get_disc();
  value.resize(num_dims);
  gradient.resize(num_dims * num_dims);
}

void Displacement<FADT>::gather(apf::MeshElement* me) {
  auto ent = apf::getMeshEntity(me);
  int num_eqs = disc->get_num_eqs();
  num_nodes = disc->get_num_nodes(ent);
  num_dofs = disc->get_num_dofs(ent);
  node_fadt.resize(num_dims * num_nodes);
  residual.resize(num_dims * num_nodes);
  elem = apf::createElement(field, me);
  apf::getVectorNodes(elem, node_st);
  for (int n = 0; n < num_nodes; ++n) {
    for (int d = 0; d < num_dims; ++d) {
      int eq = n*num_eqs + d;
      nodal(n, d).diff(eq, num_dofs);
      nodal(n, d).val() = node_st[n][d];
      resid(n, d) = 0.0;
    }
  }
}

void Displacement<FADT>::at_point(apf::Vector3 const& p, double, double) {
  auto me = apf::getMeshElement(elem);
  apf::getBF(shape, me, p, BF);
  apf::getGradBF(shape, me, p, GBF);
  for (int i = 0; i < num_dims; ++i) {
    val(i) = nodal(0, i) * BF[0];
    for (int n = 1; n < num_nodes; ++n)
      val(i) += nodal(n, i) * BF[n];
    for (int j = 0; j < num_dims; ++j) {
      grad(i, j) = nodal(0, i) * GBF[0][j];
      for (int n = 1; n < num_nodes; ++n)
        grad(i, j) += nodal(n, i) * GBF[n][j];
    }
  }
}

void Displacement<FADT>::scatter_none(SolInfo*) {
}

void Displacement<FADT>::scatter_primal(SolInfo* s) {
  using Teuchos::arrayView;
  auto ent = apf::getMeshEntity(elem);
  auto R = s->ghost->R->get1dViewNonConst();
  auto dRdu = s->ghost->dRdu;
  std::vector<LO> cols;
  disc->get_lids(ent, cols);
  auto c = arrayView(&cols[0], num_dofs);
  for (int n = 0; n < num_nodes; ++n) {
    for (int d = 0; d < num_dims; ++d) {
      auto v = resid(n, d);
      auto view = arrayView(&(v.fastAccessDx(0)), num_dofs);
      LO row = disc->get_lid(ent, n, d);
      R[row] += v.val();
      dRdu->sumIntoLocalValues(row, c, view, num_dofs);
    }
  }
}

void Displacement<FADT>::scatter_adjoint(SolInfo* s) {
  using Teuchos::arrayView;
  auto ent = apf::getMeshEntity(elem);
  auto R = s->ghost->R->get1dViewNonConst();
  auto dRduT = s->ghost->dRdu;
  std::vector<LO> cols;
  disc->get_lids(ent, cols);
  for (int n = 0; n < num_nodes; ++n) {
    for (int d = 0; d < num_dims; ++d) {
      auto v = resid(n, d);
      auto view = arrayView(&(v.fastAccessDx(0)), num_dofs);
      LO row = disc->get_lid(ent, n, d);
      R[row] += v.val();
      for (int dof = 0; dof < num_dofs; ++dof)
        dRduT->sumIntoLocalValues(
            cols[dof], arrayView(&row, 1), arrayView(&view[dof], 1));
    }
  }
}

void Displacement<FADT>::scatter(SolInfo* s) {
  op(this, s);
  apf::destroyElement(elem);
  elem = 0;
}

void Displacement<FADT>::post_process(SolInfo*) {
  disc = 0;
}

}
