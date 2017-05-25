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

void set_dbc_values(Physics* physics, const double t) {
  auto params = physics->get_dbc_params();
  auto indexer = physics->get_indexer();
  validate_params(indexer, params);
  for (auto it = params.begin(); it != params.end(); ++it) {
    auto entry = params.entry(it);
    auto a = getValue<Array<std::string> >(entry);
    auto fld = a[0];
    auto set = a[1];
    auto val = a[2];
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
    apf::synchronize(apf_f);
  }
}

template <> void apply_pre_primal_dbcs<goal::Traits::Residual>(
    Physics* physics, SolInfo* info) {
  auto params = physics->get_dbc_params();
  auto indexer = physics->get_indexer();
  auto R = info->owned->R->get1dViewNonConst();
  for (auto it = params.begin(); it != params.end(); ++it) {
    auto param_entry = params.entry(it);
    auto a = getValue<Array<std::string> >(param_entry);
    auto fld = a[0];
    auto set = a[1];
    auto idx = indexer->get_field_idx(fld);
    auto nodes = indexer->get_node_set_nodes(set);
    for (size_t node = 0; node < nodes.size(); ++node) {
      LO row = indexer->get_owned_lid(idx, nodes[node]);
      R[row] = 0.0;
    }
  }
}

template <> void apply_pre_primal_dbcs<goal::Traits::Jacobian>(
    Physics* physics, SolInfo* info) {

  auto params = physics->get_dbc_params();
  auto indexer = physics->get_indexer();
  auto dRdu = info->owned->dRdu;
  auto transposer = rcp(new Transposer(dRdu));
  auto dRduT = transposer->createTranspose();
  auto R = info->owned->R->get1dViewNonConst();

  size_t num_entries;
  Array<LO> col_indices;
  Array<ST> col_entries;

  Array<GO> row_indices;
  Array<ST> row_entries;
  Array<GO> col_index(1);
  Array<ST> entry(1);
  entry[0] = 0.0;


  for (auto it = params.begin(); it != params.end(); ++it) {

    auto param_entry = params.entry(it);
    auto a = getValue<Array<std::string> >(param_entry);
    auto fld = a[0];
    auto set = a[1];
    auto idx = indexer->get_field_idx(fld);
    auto nodes = indexer->get_node_set_nodes(set);

    for (size_t node = 0; node < nodes.size(); ++node) {

      LO row = indexer->get_owned_lid(idx, nodes[node]);
      num_entries = dRdu->getNumEntriesInLocalRow(row);
      col_indices.resize(num_entries);
      col_entries.resize(num_entries);
      dRdu->getLocalRowCopy(row, col_indices(), col_entries(), num_entries);
      for (size_t c = 0; c < num_entries; ++c)
        if (col_indices[c] != row)
          col_entries[c] = 0.0;
      dRdu->replaceLocalValues(row, col_indices(), col_entries());

      GO gid = indexer->get_owned_map()->getGlobalElement(row);
      col_index[0] = gid;
      num_entries = dRduT->getNumEntriesInGlobalRow(gid);
      row_indices.resize(num_entries);
      row_entries.resize(num_entries);
      dRduT->getGlobalRowCopy(gid, row_indices(), row_entries(), num_entries);
      for (size_t r = 0; r < num_entries; ++r)
        if (row_indices[r] != gid)
          dRdu->replaceGlobalValues(row_indices[r], col_index(), entry());

      R[row] = 0.0;
    }
  }
}

template <> void apply_post_primal_dbcs<goal::Traits::Residual>(
    Physics* physics, SolInfo* info, const double t) {
  auto params = physics->get_dbc_params();
  auto indexer = physics->get_indexer();
  auto R = info->owned->R->get1dViewNonConst();
  for (auto it = params.begin(); it != params.end(); ++it) {
    auto param_entry = params.entry(it);
    auto a = getValue<Array<std::string> >(param_entry);
    auto fld = a[0];
    auto set = a[1];
    auto val = a[2];
    auto idx = indexer->get_field_idx(fld);
    auto f = indexer->get_field(idx);
    auto apf_f = f->get_apf_field();
    auto nodes = indexer->get_node_set_nodes(set);
    for (size_t node = 0; node < nodes.size(); ++node) {
      auto n = nodes[node];
      LO row = indexer->get_owned_lid(idx, n);
      double v = get_bc_val(f, val, n, t);
      double sol = apf::getScalar(apf_f, n.entity, n.node);
      R[row] = sol - v;
    }
  }
}

template <> void apply_post_primal_dbcs<goal::Traits::Jacobian>(
    Physics* physics, SolInfo* info, const double t) {

  auto params = physics->get_dbc_params();
  auto indexer = physics->get_indexer();
  auto dRdu = info->owned->dRdu;
  auto R = info->owned->R->get1dViewNonConst();

  size_t num_entries;
  Array<LO> col_indices;
  Array<ST> col_entries;
  Array<LO> col_index(1);
  Array<ST> diag_entry(1);
  diag_entry[0] = 1.0;

  for (auto it = params.begin(); it != params.end(); ++it) {

    auto param_entry = params.entry(it);
    auto a = getValue<Array<std::string> >(param_entry);
    auto fld = a[0];
    auto set = a[1];
    auto val = a[2];
    auto idx = indexer->get_field_idx(fld);
    auto f = indexer->get_field(idx);
    auto apf_f = f->get_apf_field();
    auto nodes = indexer->get_node_set_nodes(set);

    for (size_t node = 0; node < nodes.size(); ++node) {

      auto n = nodes[node];
      LO row = indexer->get_owned_lid(idx, n);
      col_index[0] = row;
      num_entries = dRdu->getNumEntriesInLocalRow(row);
      col_indices.resize(num_entries);
      col_entries.resize(num_entries);
      dRdu->getLocalRowCopy(row, col_indices(), col_entries(), num_entries);
      for (size_t c = 0; c < num_entries; ++c)
        col_entries[c] = 0.0;
      dRdu->replaceLocalValues(row, col_indices(), col_entries());
      dRdu->replaceLocalValues(row, col_index(), diag_entry());

      double v = get_bc_val(f, val, n, t);
      double sol = apf::getScalar(apf_f, n.entity, n.node);
      R[row] = sol - v;
    }
  }
}

} // end namespace goal
