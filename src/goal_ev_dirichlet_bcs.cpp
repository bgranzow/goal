#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_control.hpp"
#include "goal_ev_dirichlet_bcs.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_solution_info.hpp"
#include "goal_workset.hpp"

namespace goal {

using Teuchos::rcp;
using Teuchos::rcpFromRef;

static void validate_params(RCP<Indexer> indexer, RCP<const ParameterList> p) {
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;
  for (auto it = p->begin(); it != p->end(); ++it) {
    auto entry = p->entry(it);
    auto a = getValue<Array<std::string> >(entry);
    assert(a.size() == 4);
    auto idx = stoi(a[0]);
    auto set = a[2];
    indexer->get_node_set_nodes(set, idx);
  }
}

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
    apf::Node const& node, const double t, bool is_adjoint) {
  double v = 0.0;
  if (is_adjoint)
    v = 0.0;
  else if (f->get_basis_type() == LAGRANGE)
    v = get_interpolatory_bc_val(f, val, node, t);
  else if (f->get_basis_type() == HIERARCHIC)
    v = get_hierarchic_bc_val(f, val, node, t);
  return v;
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
    EvalData workset, Teuchos::Array<std::string> const& a) {
  auto idx = stoi(a[0]);
  auto cmp = stoi(a[1]);
  auto set = a[2];
  auto val = a[3];
  auto u = info->owned->u;
  auto R = info->owned->R;
  auto sol = u->get1dView();
  auto res = R->get1dViewNonConst();
  auto field = indexer->get_field(idx);
  auto t = workset.t_current;
  auto nodes = indexer->get_node_set_nodes(set, idx);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    auto node = nodes[i];
    LO row = indexer->get_owned_lid(idx, node, cmp);
    double v = get_bc_val(field, val, node, t, is_adjoint);
    res[row] = sol[row] - v;
  }
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Residual, TRAITS>::evaluateFields(
    EvalData workset) {
  using Teuchos::Array;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;
  for (auto i = params->begin(); i != params->end(); ++i) {
    auto entry = params->entry(i);
    auto a = getValue<Array<std::string> >(entry);
    apply_bc(workset, a);
  }
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
    EvalData workset, Teuchos::Array<std::string> const& a) {
  auto idx = stoi(a[0]);
  auto cmp = stoi(a[1]);
  auto set = a[2];
  auto val = a[3];
  auto u = info->owned->u;
  auto R = info->owned->R;
  auto dRdu = info->owned->dRdu;
  auto sol = u->get1dView();
  auto res = R->get1dViewNonConst();
  auto field = indexer->get_field(idx);
  auto t = workset.t_current;
  size_t num_entries;
  Teuchos::Array<LO> index(1);
  Teuchos::Array<ST> entry(1);
  Teuchos::Array<LO> indices;
  Teuchos::Array<ST> entries;
  entry[0] = 1.0;
  auto nodes = indexer->get_node_set_nodes(set, idx);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    auto node = nodes[i];
    LO row = indexer->get_owned_lid(idx, node, cmp);
    double v = get_bc_val(field, val, node, t, is_adjoint);
    res[row] = sol[row] - v;
    index[0] = row;
    num_entries = dRdu->getNumEntriesInLocalRow(row);
    indices.resize(num_entries);
    entries.resize(num_entries);
    dRdu->getLocalRowCopy(row, indices(), entries(), num_entries);
    for (size_t c = 0; c < num_entries; ++c) entries[c] = 0.0;
    dRdu->replaceLocalValues(row, indices(), entries());
    dRdu->replaceLocalValues(row, index(), entry());
  }
}

template <typename TRAITS>
void DirichletBCs<goal::Traits::Jacobian, TRAITS>::evaluateFields(
    EvalData workset) {
  using Teuchos::Array;
  using Teuchos::ParameterEntry;
  using Teuchos::getValue;
  for (auto i = params->begin(); i != params->end(); ++i) {
    auto entry = params->entry(i);
    auto a = getValue<Array<std::string> >(entry);
    apply_bc(workset, a);
  }
}

template class DirichletBCs<goal::Traits::Residual, goal::Traits>;
template class DirichletBCs<goal::Traits::Jacobian, goal::Traits>;

}  // namespace goal
