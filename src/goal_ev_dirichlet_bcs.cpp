#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_control.hpp"
#include "goal_ev_dirichlet_bcs.hpp"
#include "goal_field.hpp"
#include "goal_physics.hpp"
#include "goal_indexer.hpp"
#include "goal_solution_info.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::Array;
using Teuchos::ParameterEntry;
using Teuchos::getValue;

static double get_interpolatory_bc_val(RCP<Field> f, std::string const& val,
    apf::Node const& node, const double t) {
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

static double get_hierarchic_bc_val(RCP<Field> f, std::string const& val,
    apf::Node const& node, const double t) {
  double v = 0.0;
  auto m = f->get_apf_mesh();
  if (m->getType(node.entity) == apf::Mesh::VERTEX) {
    apf::Vector3 p;
    m->getPoint(node.entity, node.node, p);
    v = eval(val, p[0], p[1], p[2], t);
  }
  return v;
}

static double get_bc_val(RCP<Field> f, std::string const& val,
    apf::Node const& node, const double t) {
  double v = 0.0;
  if (f->get_basis_type() == LAGRANGE)
    v = get_interpolatory_bc_val(f, val, node, t);
  else if (f->get_basis_type() == HIERARCHIC)
    v = get_hierarchic_bc_val(f, val, node, t);
  return v;
}

void set_dbc_values(RCP<Physics> phy, const double t) {
  auto params = phy->get_dbc_params();
  auto indexer = phy->get_indexer();
  for (auto it = params->begin(); it != params->end(); ++it) {
    auto entry = params->entry(it);
    auto a = getValue<Array<std::string> >(entry);
    auto idx = stoi(a[0]);
    auto cmp = stoi(a[1]);
    auto set = a[2];
    auto val = a[3];
    auto f = indexer->get_field(idx);
    auto apf_f = f->get_apf_field();
    auto nodes = indexer->get_node_set_nodes(set, idx);
    double vals[3] = {0, 0, 0};
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      auto node = nodes[i];
      auto v = get_bc_val(f, val, node, t);
      auto e = node.entity;
      auto n = node.node;
      apf::getComponents(apf_f, e, n, vals);
      vals[cmp] = v;
      apf::setComponents(apf_f, e, n, vals);
    }
  }
}

static void validate_params(RCP<Indexer> indexer, RCP<const ParameterList> p) {
  for (auto it = p->begin(); it != p->end(); ++it) {
    auto entry = p->entry(it);
    auto a = getValue<Array<std::string> >(entry);
    assert(a.size() == 4);
    auto idx = stoi(a[0]);
    auto set = a[2];
    indexer->get_node_set_nodes(set, idx);
  }
}

template <typename TRAITS>
DirichletBCs<goal::Traits::Residual, TRAITS>::DirichletBCs(
    RCP<const ParameterList> p, RCP<Indexer> i, bool adj) {
  params = p;
  indexer = i;
  is_adjoint = adj;
  validate_params(indexer, params);
  auto name = "Dirichlet BCs";
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Residual, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  (void)d;
  (void)fm;
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Residual, TRAITS>::preEvaluate(PreEvalData i) {
  info = rcpFromRef(i);
  assert(info->owned->R != Teuchos::null);
  assert(info->owned->u != Teuchos::null);
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Residual, TRAITS>::apply_bc(
    Teuchos::Array<std::string> const& a) {
  auto idx = stoi(a[0]);
  auto cmp = stoi(a[1]);
  auto set = a[2];
  auto R = info->owned->R->get1dViewNonConst();
  auto nodes = indexer->get_node_set_nodes(set, idx);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    auto node = nodes[i];
    LO row = indexer->get_owned_lid(idx, node, cmp);
    R[row] = 0.0;
  }
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Residual, TRAITS>::evaluateFields(
    EvalData workset) {
  for (auto i = params->begin(); i != params->end(); ++i) {
    auto entry = params->entry(i);
    auto a = getValue<Array<std::string> >(entry);
    apply_bc(a);
  }
  (void)workset;
}

template <typename TRAITS>
DirichletBCs<goal::Traits::Jacobian, TRAITS>::DirichletBCs(
    RCP<const ParameterList> p, RCP<Indexer> i, bool adj) {
  params = p;
  indexer = i;
  is_adjoint = adj;
  validate_params(indexer, params);
  auto name = "Dirichlet BCs";
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Jacobian, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  (void)d;
  (void)fm;
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Jacobian, TRAITS>::preEvaluate(PreEvalData i) {
  info = rcpFromRef(i);
  assert(info->owned->R != Teuchos::null);
  assert(info->owned->u != Teuchos::null);
  assert(info->owned->dRdu != Teuchos::null);
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Jacobian, TRAITS>::apply_bc(
    Teuchos::Array<std::string> const& a) {
  auto idx = stoi(a[0]);
  auto cmp = stoi(a[1]);
  auto set = a[2];
  auto R = info->owned->R->get1dViewNonConst();
  auto dRdu = info->owned->dRdu;
  size_t num_entries;
  Teuchos::Array<LO> indices;
  Teuchos::Array<ST> entries;
  auto nodes = indexer->get_node_set_nodes(set, idx);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    auto node = nodes[i];
    LO row = indexer->get_owned_lid(idx, node, cmp);
    num_entries = dRdu->getNumEntriesInLocalRow(row);
    indices.resize(num_entries);
    entries.resize(num_entries);
    dRdu->getLocalRowCopy(row, indices(), entries(), num_entries);
    for (size_t c = 0; c < num_entries; ++c)
      if (indices[c] != row)
        entries[c] = 0.0;
    dRdu->replaceLocalValues(row, indices(), entries());
    R[row] = 0.0;
  }
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Jacobian, TRAITS>::evaluateFields(
    EvalData workset) {
  for (auto i = params->begin(); i != params->end(); ++i) {
    auto param_entry = params->entry(i);
    auto a = getValue<Array<std::string> >(param_entry);
    apply_bc(a);
  }
  (void)workset;
}

template class DirichletBCs<goal::Traits::Residual, goal::Traits>;
template class DirichletBCs<goal::Traits::Jacobian, goal::Traits>;

}  // namespace goal
