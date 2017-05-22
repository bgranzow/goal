#include "goal_assembly.hpp"
#include "goal_control.hpp"
#include "goal_dbcs.hpp"
#include "goal_discretization.hpp"
#include "goal_physics.hpp"
#include "goal_sol_info.hpp"
#include "goal_workset.hpp"

namespace goal {

using Residual = goal::Traits::Residual;
using Jacobian = goal::Traits::Jacobian;

static void load_elem_info(
    Workset& ws, Discretization* d, const int es_idx, const int ws_idx) {
  auto set = d->get_elem_set_name(es_idx);
  ws.set = set;
  ws.entities = d->get_elems(set, ws_idx);
  ws.size = ws.entities.size();
}

static void load_side_info(
    Workset& ws, Discretization* d, const int ss_idx, const int ws_idx) {
  auto set = d->get_side_set_name(ss_idx);
  ws.set = set;
  ws.entities = d->get_sides(set, ws_idx);
  ws.size = ws.entities.size();
}

template <typename EvalT>
void assemble_volumetric(
    Workset& ws, Physics* p, SolInfo* i, Discretization* d) {
  auto f = p->get_volumetric();
  for (int es_idx = 0; es_idx < f.size(); ++es_idx) {
    f[es_idx]->template preEvaluate<EvalT>(i);
    auto num_ws = d->get_num_elem_worksets(es_idx);
    for (int ws_idx = 0; ws_idx < num_ws; ++ws_idx) {
      load_elem_info(ws, d, es_idx, ws_idx);
      f[es_idx]->template evaluateFields<EvalT>(ws);
    }
    f[es_idx]->template postEvaluate<EvalT>(i);
  }
}

template <typename EvalT>
void assemble_neumann(
    Workset& ws, Physics* p, SolInfo* i, Discretization* d) {
  auto f = p->get_neumann();
  for (int ss_idx = 0; ss_idx < f.size(); ++ss_idx) {
    f[ss_idx]->template preEvaluate<EvalT>(i);
    auto num_ws = d->get_num_side_worksets(ss_idx);
    for (int ws_idx = 0; ws_idx < num_ws; ++ws_idx) {
      load_side_info(ws, d, ss_idx, ws_idx);
      f[ss_idx]->template evaluateFields<EvalT>(ws);
    }
    f[ss_idx]->template postEvaluate<EvalT>(i);
  }
}

template <typename EvalT>
void assemble_dirichlet(
    Workset& ws, Physics* p, SolInfo* i, Discretization* d) {
  auto f = p->get_dirichlet();
  f->template preEvaluate<EvalT>(i);
  f->template evaluateFields<EvalT>(ws);
  f->template postEvaluate<EvalT>(i);
  (void)d;
}

void compute_primal_residual(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old) {
  auto t0 = time();
  Workset ws;
  ws.t_now = t_now;
  ws.t_old = t_old;
  i->owned->R->putScalar(0.0);
  i->ghost->R->putScalar(0.0);
  assemble_volumetric<Residual>(ws, p, i, d);
  assemble_neumann<Residual>(ws, p, i, d);
  i->gather_R();
  assemble_dirichlet<Residual>(ws, p, i, d);
  apply_primal_dbcs<Residual>(p, i);
  auto t1 = time();
  print(" > residual computed in %f seconds", t1 - t0);
}

void compute_primal_jacobian(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old) {
  auto t0 = time();
  Workset ws;
  ws.t_now = t_now;
  ws.t_old = t_old;
  i->owned->R->putScalar(0.0);
  i->ghost->R->putScalar(0.0);
  i->owned->dRdu->resumeFill();
  i->ghost->dRdu->resumeFill();
  i->owned->dRdu->setAllToScalar(0.0);
  i->ghost->dRdu->setAllToScalar(0.0);
  assemble_volumetric<Jacobian>(ws, p, i, d);
  assemble_neumann<Jacobian>(ws, p, i, d);
  i->ghost->dRdu->fillComplete();
  i->gather_R();
  i->gather_dRdu();
  assemble_dirichlet<Jacobian>(ws, p, i, d);
  apply_primal_dbcs<Jacobian>(p, i, true);
  i->owned->dRdu->fillComplete();
  auto t1 = time();
  print(" > jacobian computed in %f seconds", t1 - t0);
}

void compute_dual_jacobian(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old) {
  auto t0 = time();
  Workset ws;
  ws.t_now = t_now;
  ws.t_old = t_old;
  i->owned->R->putScalar(0.0);
  i->owned->R->putScalar(0.0);
  i->owned->dJdu->putScalar(0.0);
  i->ghost->dJdu->putScalar(0.0);
  i->owned->dRdu->resumeFill();
  i->ghost->dRdu->resumeFill();
  i->owned->dRdu->setAllToScalar(0.0);
  i->ghost->dRdu->setAllToScalar(0.0);
  assemble_volumetric<Jacobian>(ws, p, i, d);
  assemble_neumann<Jacobian>(ws, p, i, d);
  i->ghost->dRdu->fillComplete();
  i->gather_R();
  i->gather_dJdu();
  i->gather_dRdu();
  assemble_dirichlet<Jacobian>(ws, p, i, d);
  apply_dual_dbcs<Jacobian>(p, i, true);
  auto t1 = time();
  print(" > jacobian computed in %f seconds", t1 - t0);
}

void compute_error_residual(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old) {
  auto t0 = time();
  Workset ws;
  ws.t_now = t_now;
  ws.t_old = t_old;
  i->owned->R->putScalar(0.0);
  i->ghost->R->putScalar(0.0);
  assemble_volumetric<Residual>(ws, p, i, d);
  assemble_neumann<Residual>(ws, p, i, d);
  i->gather_R();
  assemble_dirichlet<Residual>(ws, p, i, d);
  apply_primal_dbcs<Residual>(p, i, t_now);
  auto t1 = time();
  print(" > residual computed in %f seconds", t1 - t0);
}

template void assemble_volumetric<Residual>(Workset& ws, Physics* p, SolInfo* i, Discretization* d);
template void assemble_volumetric<Jacobian>(Workset& ws, Physics* p, SolInfo* i, Discretization* d);
template void assemble_neumann<Residual>(Workset& ws, Physics* p, SolInfo* i, Discretization* d);
template void assemble_neumann<Jacobian>(Workset& ws, Physics* p, SolInfo* i, Discretization* d);
template void assemble_dirichlet<Residual>(Workset& ws, Physics* p, SolInfo* i, Discretization* d);
template void assemble_dirichlet<Jacobian>(Workset& ws, Physics* p, SolInfo* i, Discretization* d);

} // end namespace goal
