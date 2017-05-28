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

static apf::FieldShape* get_basis(int basis, int p) {
  apf::FieldShape* shape = 0;
  if (basis == LAGRANGE) {
    if (p == 2)
      shape = apf::getSerendipity();
    else
      shape = apf::getLagrange(p);
  } else if (basis == HIERARCHICAL)
    shape = apf::getHierarchic(p);
  return shape;
}

Field::Field(FieldInfo const& info) {
  disc = info.disc;
  myname = info.name;
  p_order = info.p_order;
  q_degree = info.q_degree;
  basis_type = info.basis_type;
  apf_mesh = disc->get_apf_mesh();
  apf_basis = get_basis(basis_type, p_order);
  auto n = myname.c_str();
  auto t = apf::SCALAR;
  apf_field = apf::createField(apf_mesh, n, t, apf_basis);
  apf::zeroField(apf_field);
  idx = -1;
  seed = 1.0;
}

Field::~Field() {
  if (apf_field)
    apf::destroyField(apf_field);
}

int Field::get_num_dims() const {
  return disc->get_num_dims();
}

int Field::get_num_nodes(const int t) const {
  auto es = apf_basis->getEntityShape(t);
  return es->countNodes();
}

int Field::get_num_vtx(const int t) const {
  return apf::Mesh::adjacentCount[t][0];
}

int Field::get_num_ips(const int t) const {
  return apf::countGaussPoints(t, q_degree);
}

std::string Field::name() const {
  return myname;
}

std::string Field::g_name() const {
  return "grad_" + myname;
}

std::string Field::resid_name() const {
  return myname + "_resid";
}

std::string Field::basis_name() const {
  std::map<int, std::string> map = {
    {LAGRANGE, "lagrange"}, {HIERARCHICAL, "hierarchical"}};
  std::ostringstream oss;
  oss << map[basis_type] << p_order;
  auto n = oss.str();
  return n;
}

std::string Field::g_basis_name() const {
  return "grad_" + basis_name();
}

std::string Field::wdv_name() const {
  return basis_name() + "_wdv";
}

RCP<PHX::DataLayout> Field::dl(const int t) {
  int ws = disc->get_ws_size();
  int n = get_num_nodes(t);
  return rcp(new PHX::MDALayout<Ent, Node>(ws, n));
}

RCP<PHX::DataLayout> Field::g_dl(const int t) {
  int ws = disc->get_ws_size();
  int n = get_num_nodes(t);
  int d = get_num_dims();
  return rcp(new PHX::MDALayout<Ent, Node, Dim>(ws, n, d));
}

RCP<PHX::DataLayout> Field::w_dl(const int t) {
  int ws = disc->get_ws_size();
  int n = get_num_nodes(t);
  int ip = get_num_ips(t);
  return rcp(new PHX::MDALayout<Ent, Node, IP>(ws, n, ip));
}

RCP<PHX::DataLayout> Field::g_w_dl(const int t) {
  int ws = disc->get_ws_size();
  int n = get_num_nodes(t);
  int ip = get_num_ips(t);
  int d = get_num_dims();
  return rcp(new PHX::MDALayout<Ent, Node, IP, Dim>(ws, n, ip, d));
}

RCP<PHX::DataLayout> Field::ip_dl(const int t) {
  int ws = disc->get_ws_size();
  int ip = get_num_ips(t);
  return rcp(new PHX::MDALayout<Ent, IP>(ws, ip));
}

RCP<PHX::DataLayout> Field::g_ip_dl(const int t) {
  int ws = disc->get_ws_size();
  int ip = get_num_ips(t);
  int d = get_num_dims();
  return rcp(new PHX::MDALayout<Ent, IP, Dim>(ws, ip, d));
}

RCP<PHX::DataLayout> Field::PU_dl(const int t) {
  int ws = disc->get_ws_size();
  int v = get_num_vtx(t);
  int ip = get_num_ips(t);
  return rcp(new PHX::MDALayout<Ent, Node, IP>(ws, v, ip));
}

RCP<PHX::DataLayout> Field::g_PU_dl(const int t) {
  int ws = disc->get_ws_size();
  int v = get_num_vtx(t);
  int ip = get_num_ips(t);
  int d = get_num_dims();
  return rcp(new PHX::MDALayout<Ent, Node, IP, Dim>(ws, v, ip, d));
}

RCP<PHX::DataLayout> Field::ip0_dl(const int t) {
  int ws = disc->get_ws_size();
  int ip = get_num_ips(t);
  return rcp(new PHX::MDALayout<Ent, IP>(ws, ip));
}

RCP<PHX::DataLayout> Field::ip1_dl(const int t) {
  int ws = disc->get_ws_size();
  int ip = get_num_ips(t);
  int d = get_num_dims();
  return rcp(new PHX::MDALayout<Ent, IP, Dim>(ws, ip, d));
}

RCP<PHX::DataLayout> Field::ip2_dl(const int t) {
  int ws = disc->get_ws_size();
  int ip = get_num_ips(t);
  int d = get_num_dims();
  return rcp(new PHX::MDALayout<Ent, IP, Dim, Dim>(ws, ip, d, d));
}

RCP<PHX::DataLayout> Field::ent0_dl(const int t) {
  int ws = disc->get_ws_size();
  (void)t;
  return rcp(new PHX::MDALayout<Ent>(ws));
}

Field* create_field(FieldInfo const& i) {
  return new Field(i);
}

void destroy_field(Field* f) {
  delete f;
}

} // end namespace goal
