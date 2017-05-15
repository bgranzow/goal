#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_field.hpp>
#include <Phalanx_DataLayout_MDALayout.hpp>
#include <Teuchos_ParameterList.hpp>

namespace test {

static const int ws = 1000;
static const std::string name = "u";
static const int p = 2;
static const int q = 2;
static const int idx = 0;
static const int basis_type = goal::LAGRANGE;
static const double seed = 1.0;

static goal::Discretization* load_disc(char** argv) {
  Teuchos::ParameterList p;
  p.set<std::string>("geom file", argv[1]);
  p.set<std::string>("mesh file", argv[2]);
  p.set<std::string>("assoc file", argv[3]);
  p.set<bool>("reorder mesh", true);
  p.set<int>("workset size", ws);
  return goal::create_disc(p);
}

static goal::Field* make_field(goal::Discretization* d) {
  goal::FieldInfo i = { d, name, p, q, basis_type };
  auto f = goal::create_field(i);
  f->set_seed(seed);
  f->set_associated_dof_idx(idx);
  return f;
}

static void check_dl(goal::Field* f, const int t) {
  auto dl = f->dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 2);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_nodes(t));
}

static void check_g_dl(goal::Field* f, const int t) {
  auto dl = f->g_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 3);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_nodes(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(2) == f->get_num_dims());
}

static void check_w_dl(goal::Field* f, const int t) {
  auto dl = f->w_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 3);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_nodes(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(2) == f->get_num_ips(t));
}

static void check_g_w_dl(goal::Field* f, const int t) {
  auto dl = f->g_w_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 4);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_nodes(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(2) == f->get_num_ips(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(3) == f->get_num_dims());
}

static void check_ip_dl(goal::Field* f, const int t) {
  auto dl = f->ip_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 2);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_ips(t));
}

static void check_g_ip_dl(goal::Field* f, const int t) {
  auto dl = f->g_ip_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 3);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_ips(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(2) == f->get_num_dims());
}

static void check_PU_dl(goal::Field* f, const int t) {
  auto dl = f->PU_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 3);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_vtx(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(2) == f->get_num_ips(t));
}

static void check_g_PU_dl(goal::Field* f, const int t) {
  auto dl = f->g_PU_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 4);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_vtx(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(2) == f->get_num_ips(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(3) == f->get_num_dims());
}

static void check_ip0_dl(goal::Field* f, const int t) {
  auto dl = f->ip0_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 2);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_ips(t));
}

static void check_ip1_dl(goal::Field* f, const int t) {
  auto dl = f->ip1_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 3);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_ips(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(2) == f->get_num_dims());
}

static void check_ip2_dl(goal::Field* f, const int t) {
  auto dl = f->ip2_dl(t);
  GOAL_ALWAYS_ASSERT((int)dl->rank() == 4);
  GOAL_ALWAYS_ASSERT(dl->dimension(0) == ws);
  GOAL_ALWAYS_ASSERT((int)dl->dimension(1) == f->get_num_ips(t));
  GOAL_ALWAYS_ASSERT((int)dl->dimension(2) == f->get_num_dims());
  GOAL_ALWAYS_ASSERT((int)dl->dimension(3) == f->get_num_dims());
}

static void check_all_dl(goal::Field* f, const int t) {
  goal::print("checking dl's for type: %d", t);
  check_dl(f, t);
  check_g_dl(f, t);
  check_w_dl(f, t);
  check_g_w_dl(f, t);
  check_ip_dl(f, t);
  check_g_ip_dl(f, t);
  check_PU_dl(f, t);
  check_g_PU_dl(f, t);
  check_ip0_dl(f, t);
  check_ip1_dl(f, t);
  check_ip2_dl(f, t);
}

static void check_field(goal::Discretization* d, goal::Field* f) {
  GOAL_ALWAYS_ASSERT(p == f->get_p_order());
  GOAL_ALWAYS_ASSERT(q == f->get_q_degree());
  GOAL_ALWAYS_ASSERT(basis_type == f->get_basis_type());
  GOAL_ALWAYS_ASSERT(seed == f->get_seed_value());
  GOAL_ALWAYS_ASSERT(idx == f->get_associated_dof_idx());
  GOAL_ALWAYS_ASSERT(d->get_num_dims() == f->get_num_dims());
  goal::print("name:          %s", (f->name()).c_str());
  goal::print("g_name:        %s", (f->g_name()).c_str());
  goal::print("resid_name:    %s", (f->resid_name()).c_str());
  goal::print("basis_name:    %s", (f->basis_name()).c_str());
  goal::print("g_basis_name:  %s", (f->g_basis_name()).c_str());
  goal::print("wdv_name:      %s", (f->wdv_name()).c_str());
  for (int i = 0; i < d->get_num_elem_sets(); ++i) {
    for (int ws = 0; ws < d->get_num_elem_worksets(i); ++ws)
      check_all_dl(f, d->get_elem_type(i));
  }
  for (int i = 0; i < d->get_num_side_sets(); ++i) {
    for (int ws = 0; ws < d->get_num_side_worksets(i); ++ws)
      check_all_dl(f, d->get_side_type(i));
  }
}

} // end namespace test

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: field");
  GOAL_ALWAYS_ASSERT(argc == 4);
  auto d = test::load_disc(argv);
  auto f = test::make_field(d);
  test::check_field(d, f);
  goal::destroy_field(f);
  goal::destroy_disc(d);
  goal::finalize();
}
