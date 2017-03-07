#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <Phalanx_DataLayout_MDALayout.hpp>

#include "goal_control.hpp"
#include "goal_dimension.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"

namespace goal {

using Teuchos::rcp;

apf::FieldShape* get_basis(int basis, int p) {
  apf::FieldShape* shape = 0;
  if (basis == LAGRANGE) {
    if (p == 2)
      shape = apf::getSerendipity();
    else
      shape = apf::getLagrange(p);
  } else if (basis == HIERARCHIC)
    shape = apf::getHierarchic(p);
  return shape;
}

Field::Field(FieldInfo* info) {
  disc = info->disc;
  name = info->name;
  p_order = info->p_order;
  q_degree = info->q_degree;
  value_type = info->value_type;
  basis_type = info->basis_type;
  apf_mesh = disc->get_apf_mesh();
  apf_basis = get_basis(basis_type, p_order);
  auto n = name.c_str();
  apf_field = apf::createField(apf_mesh, n, value_type, apf_basis);
  apf::zeroField(apf_field);
  idx = -1;
  seed = 1.0;
  set_elem_block(0);
}

Field::~Field() {
  if (apf_field) apf::destroyField(apf_field);
}

void Field::set_elem_block(const int idx) {
  elem_type = disc->get_elem_type(idx);
}

int Field::get_num_dims() const { return disc->get_num_dims(); }

int Field::get_num_components() const {
  int nc = 0;
  if (value_type == SCALAR)
    nc = 1;
  else if (value_type == VECTOR)
    nc = disc->get_num_dims();
  else
    fail("unsupported value type");
  return nc;
}

int Field::get_num_elem_nodes() const {
  auto es = apf_basis->getEntityShape(elem_type);
  return es->countNodes();
}

int Field::get_num_elem_vtx() const {
  return apf::Mesh::adjacentCount[elem_type][0];
}

int Field::get_num_elem_dofs() const {
  auto nc = get_num_components();
  auto nn = get_num_elem_nodes();
  return nc * nn;
}

int Field::get_num_elem_ips() const {
  return apf::countGaussPoints(elem_type, q_degree);
}

std::string Field::get_name() const { return name; }

std::string Field::get_grad_name() const { return "grad_" + name; }

std::string Field::get_basis_name() const {
  std::map<int, std::string> map = {
      {LAGRANGE, "lagrange"}, {HIERARCHIC, "hierarchic"}};
  std::ostringstream oss;
  oss << map[get_value_type()] << get_p_order();
  auto name = oss.str();
  return name;
}

std::string Field::get_grad_basis_name() const {
  return "grad_" + get_basis_name();
}

std::string Field::get_residual_name() const {
  assert(dof);
  return name + "_residual";
}

std::string Field::get_wdv_name() const { return name + "_wdv"; }

RCP<PHX::DataLayout> Field::get_dl() {
  int ws_size = disc->get_ws_size();
  int n_nodes = get_num_elem_nodes();
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(new PHX::MDALayout<Elem, Node>(ws_size, n_nodes));
  else if (value_type == VECTOR)
    return rcp(new PHX::MDALayout<Elem, Node, Dim>(ws_size, n_nodes, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_grad_dl() {
  int ws_size = disc->get_ws_size();
  int n_nodes = get_num_elem_nodes();
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(new PHX::MDALayout<Elem, Node, Dim>(ws_size, n_nodes, n_dims));
  else if (value_type == VECTOR)
    return rcp(new PHX::MDALayout<Elem, Node, Dim, Dim>(
        ws_size, n_nodes, n_dims, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_weight_dl() {
  int ws_size = disc->get_ws_size();
  int n_nodes = get_num_elem_nodes();
  int n_ips = get_num_elem_ips();
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(new PHX::MDALayout<Elem, Node, IP>(ws_size, n_nodes, n_ips));
  else if (value_type == VECTOR)
    return rcp(new PHX::MDALayout<Elem, Node, IP, Dim>(
        ws_size, n_nodes, n_ips, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_grad_weight_dl() {
  int ws_size = disc->get_ws_size();
  int n_nodes = get_num_elem_nodes();
  int n_ips = get_num_elem_ips();
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(new PHX::MDALayout<Elem, Node, IP, Dim>(
        ws_size, n_nodes, n_ips, n_dims));
  else if (value_type == VECTOR)
    return rcp(new PHX::MDALayout<Elem, Node, IP, Dim, Dim>(
        ws_size, n_nodes, n_ips, n_dims, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_interpolated_dl() {
  int ws_size = disc->get_ws_size();
  int n_ips = get_num_elem_ips();
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(new PHX::MDALayout<Elem, IP>(ws_size, n_ips));
  else if (value_type == VECTOR)
    return rcp(new PHX::MDALayout<Elem, IP, Dim>(ws_size, n_ips, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_grad_interpolated_dl() {
  int ws_size = disc->get_ws_size();
  int n_ips = get_num_elem_ips();
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(new PHX::MDALayout<Elem, IP, Dim>(ws_size, n_ips, n_dims));
  else if (value_type == VECTOR)
    return rcp(
        new PHX::MDALayout<Elem, IP, Dim, Dim>(ws_size, n_ips, n_dims, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_PU_dl() {
  int ws_size = disc->get_ws_size();
  int n_vtx = apf::Mesh::adjacentCount[elem_type][0];
  int n_ips = get_num_elem_ips();
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(new PHX::MDALayout<Elem, Node, IP>(ws_size, n_vtx, n_ips));
  else if (value_type == VECTOR)
    return rcp(
        new PHX::MDALayout<Elem, Node, IP, Dim>(ws_size, n_vtx, n_ips, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_grad_PU_dl() {
  int ws_size = disc->get_ws_size();
  int n_vtx = apf::Mesh::adjacentCount[elem_type][0];
  int n_ips = get_num_elem_ips();
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(
        new PHX::MDALayout<Elem, Node, IP, Dim>(ws_size, n_vtx, n_ips, n_dims));
  else if (value_type == VECTOR)
    return rcp(new PHX::MDALayout<Elem, Node, IP, Dim, Dim>(
        ws_size, n_vtx, n_ips, n_dims, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_residual_PU_dl() {
  int ws_size = disc->get_ws_size();
  int n_vtx = apf::Mesh::adjacentCount[elem_type][0];
  int n_dims = get_num_dims();
  if (value_type == SCALAR)
    return rcp(new PHX::MDALayout<Elem, Node>(ws_size, n_vtx));
  else if (value_type == VECTOR)
    return rcp(new PHX::MDALayout<Elem, Node, Dim>(ws_size, n_vtx, n_dims));
  else
    return Teuchos::null;
}

RCP<PHX::DataLayout> Field::get_scalar_ip_dl() {
  int ws_size = disc->get_ws_size();
  int n_ips = get_num_elem_ips();
  return rcp(new PHX::MDALayout<Elem, IP>(ws_size, n_ips));
}

RCP<PHX::DataLayout> Field::get_vector_ip_dl() {
  int ws_size = disc->get_ws_size();
  int n_ips = get_num_elem_ips();
  int n_dims = get_num_dims();
  return rcp(new PHX::MDALayout<Elem, IP, Dim>(ws_size, n_ips, n_dims));
}

RCP<PHX::DataLayout> Field::get_tensor_ip_dl() {
  int ws_size = disc->get_ws_size();
  int n_ips = get_num_elem_ips();
  int n_dims = get_num_dims();
  return rcp(
      new PHX::MDALayout<Elem, IP, Dim, Dim>(ws_size, n_ips, n_dims, n_dims));
}

RCP<PHX::DataLayout> Field::get_elem_scalar_dl() {
  int ws_size = disc->get_ws_size();
  return rcp(new PHX::MDALayout<Elem>(ws_size));
}

void project_field(RCP<Field> to, RCP<Field> from) {
  auto to_value_type = to->get_value_type();
  auto from_value_type = from->get_value_type();
  auto to_basis_type = to->get_basis_type();
  auto from_basis_type = from->get_basis_type();
  assert(to_value_type == from_value_type);
  assert(to_basis_type == from_basis_type);
  auto apf_to = to->get_apf_field();
  auto apf_from = from->get_apf_field();
  if (to_basis_type == LAGRANGE)
    apf::projectField(apf_to, apf_from);
  else if (to_basis_type == HIERARCHIC)
    apf::projectHierarchicField(apf_to, apf_from);
  else
    fail("project: unknown value type");
}

}  // namespace goal
