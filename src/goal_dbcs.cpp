#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <Tpetra_RowMatrixTransposer_decl.hpp>

#include "goal_control.hpp"
#include "goal_dbcs.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_physics.hpp"
#include "goal_sol_info.hpp"

namespace goal {

using Teuchos::rcp;
using Teuchos::Array;
using Teuchos::ParameterEntry;
using Teuchos::getValue;
using Transposer = Tpetra::RowMatrixTransposer<ST, LO, GO, KNode>;

static void validate_params(Indexer* indexer, ParameterList const& p) {
  for (auto it = p.begin(); it != p.end(); ++it) {
    auto entry = p.entry(it);
    auto a = getValue<Array<std::string> >(entry);
    GOAL_ALWAYS_ASSERT(a.size() == 3);
    auto set = a[1];
    auto nodes = indexer->get_node_set_nodes(set);
  }
}

static double get_interp_val(
    Field* f,
    std::string const& val,
    apf::Node const& node,
    const double t) {
  apf::Vector3 p;
  apf::Vector3 xi;
  auto e = node.entity;
  auto m = f->get_apf_mesh();
  auto type = m->getType(e);
  auto me = apf::createMeshElement(m, e);
  f->get_apf_basis()->getNodeXi(type, node.node, xi);
  apf::mapLocalToGlobal(me, xi, p);
  apf::destroyMeshElement(me);
  double v = eval(val, p[0], p[1], p[2], t);
  return v;
}

static double get_hierarchic_val(
    Field* f,
    std::string const& val,
    apf::Node const& node,
    const double t) {
  double v = 0.0;
  auto m = f->get_apf_mesh();
  if (m->getType(node.entity) == apf::Mesh::VERTEX) {
    apf::Vector3 p;
    m->getPoint(node.entity, node.node, p);
    v = eval(val, p[0], p[1], p[2], t);
  }
  return v;
}

static double get_bc_val(
    Field* f,
    std::string const& val,
    apf::Node const& node,
    const double t) {
  double v = 0.0;
  if (f->get_basis_type() == LAGRANGE)
    v = get_interp_val(f, val, node, t);
  else if (f->get_basis_type() == HIERARCHICAL)
    v = get_hierarchic_val(f, val, node, t);
  return v;
}

enum Mode { RES, JAC, QOI };

static void apply_bc_res(
    int idx,
    SolInfo* info,
    Indexer* indexer,
    std::vector<apf::Node> const& nodes) {
  auto R = info->owned->R;
  for (size_t n = 0; n < nodes.size(); ++n) {
    LO row = indexer->get_owned_lid(idx, nodes[n]);
    R->replaceLocalValue(row, 0.0);
  }
}

static void apply_bc_jac(
    int idx,
    SolInfo* info,
    Indexer* indexer,
    std::vector<apf::Node> const& nodes) {
  auto R = info->owned->R;
  auto dRdu = info->owned->dRdu;
  size_t num_entries;
  Array<LO> indices;
  Array<ST> entries;
  for (size_t n = 0; n < nodes.size(); ++n) {
    LO row = indexer->get_owned_lid(idx, nodes[n]);
    num_entries = dRdu->getNumEntriesInLocalRow(row);
    indices.resize(num_entries);
    entries.resize(num_entries);
    dRdu->getLocalRowCopy(row, indices(), entries(), num_entries);
    for (size_t c = 0; c < num_entries; ++c)
      if (indices[c] != row)
        entries[c] = 0.0;
    dRdu->replaceLocalValues(row, indices(), entries());
    R->replaceLocalValue(row, 0.0);
  }
}

static void apply_bc_qoi(
    int idx,
    SolInfo* info,
    Indexer* indexer,
    std::vector<apf::Node> const& nodes) {
  auto dJdu = info->owned->dJdu;
  for (size_t n = 0; n < nodes.size(); ++n) {
    LO row = indexer->get_owned_lid(idx, nodes[n]);
    for (size_t c = 0; c < dJdu->getNumVectors(); ++c)
      dJdu->replaceLocalValue(row, c, 0.0);
  }
}

static void apply_bc(
    SolInfo* info, Indexer* indexer, ParameterList const& p, int m) {
  for (auto it = p.begin(); it != p.end(); ++it) {
    auto entry = p.entry(it);
    auto a = getValue<Array<std::string> >(entry);
    auto fld = a[0];
    auto set = a[1];
    auto idx = indexer->get_field_idx(fld);
    auto nodes = indexer->get_node_set_nodes(set);
    switch (m) {
      case RES: apply_bc_res(idx, info, indexer, nodes); break;
      case JAC: apply_bc_jac(idx, info, indexer, nodes); break;
      case QOI: apply_bc_qoi(idx, info, indexer, nodes); break;
    }
  }
}

static void condense(
    SolInfo* info, Indexer* indexer, ParameterList const& p) {
  auto dRdu = info->owned->dRdu;
  auto transposer = rcp(new Transposer(dRdu));
  auto dRduT = transposer->createTranspose();
  size_t num_entries;
  Array<GO> indices;
  Array<ST> entries;
  Array<GO> index(1);
  Array<ST> entry(1);
  entry[0] = 0.0;
  for (auto i = p.begin(); i != p.end(); ++i) {
    auto param_entry = p.entry(i);
    auto a = getValue<Array<std::string> >(param_entry);
    auto fld = a[0];
    auto set = a[1];
    auto idx = indexer->get_field_idx(fld);
    auto nodes = indexer->get_node_set_nodes(set);
    for (std::size_t n = 0; n < nodes.size(); ++n) {
      LO dof_lid = indexer->get_owned_lid(idx, nodes[n]);
      GO dof_gid = indexer->get_owned_map()->getGlobalElement(dof_lid);
      index[0] = dof_gid;
      num_entries = dRduT->getNumEntriesInGlobalRow(dof_gid);
      indices.resize(num_entries);
      entries.resize(num_entries);
      dRduT->getGlobalRowCopy(dof_gid, indices(), entries(), num_entries);
      for (size_t r = 0; r < num_entries; ++r)
        if (indices[r] != dof_gid)
          dRdu->replaceGlobalValues(indices[r], index(), entry());
    }
  }
}

void set_dbc_values(Physics* p, const double t) {
  auto params = p->get_dbc_params();
  auto indexer = p->get_indexer();
  validate_params(indexer, params);
  for (auto it = params.begin(); it != params.end(); ++it) {
    auto entry = params.entry(it);
    auto a = getValue<Array<std::string> >(entry);
    auto fld = a[0];
    auto set = a[1];
    auto val = a[3];
    auto idx = indexer->get_field_idx(fld);
    auto f = indexer->get_field(idx);
    auto apf_f = f->get_apf_field();
    auto nodes = indexer->get_node_set_nodes(set);
    for (size_t i =0; i < nodes.size(); ++i) {
      auto node = nodes[i];
      auto v = get_bc_val(f, val, node, t);
      auto e = node.entity;
      auto n = node.node;
      apf::setScalar(apf_f, e, n, v);
    }
  }
}

template <> void apply_primal_dbcs<goal::Traits::Residual>(
    Physics* p, SolInfo* i, bool c) {
  auto params = p->get_dbc_params();
  auto indexer = p->get_indexer();
  apply_bc(i, indexer, params, RES);
  (void)c;
}

template <> void apply_primal_dbcs<goal::Traits::Jacobian>(
    Physics* p, SolInfo* i, bool c) {
  auto params = p->get_dbc_params();
  auto indexer = p->get_indexer();
  apply_bc(i, indexer, params, JAC);
  if (c) condense(i, indexer, params);
}

template <> void apply_dual_dbcs<goal::Traits::Jacobian>(
    Physics* p, SolInfo* i, bool c) {
  auto params = p->get_dbc_params();
  auto indexer = p->get_indexer();
  apply_bc(i, indexer, params, JAC);
  apply_bc(i, indexer, params, QOI);
  if (c) condense(i, indexer, params);
}

} // end namespace goal
