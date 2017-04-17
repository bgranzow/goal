#include "goal_assembly.hpp"
#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_physics.hpp"
#include "goal_solution_info.hpp"
#include "goal_workset.hpp"

namespace goal {

static void load_disc_info(
    Workset& ws, RCP<Discretization> d, const int block_idx, const int ws_idx) {
  auto block = d->get_elem_block_name(block_idx);
  ws.block = block;
  ws.entities = d->get_elems(block, ws_idx);
  ws.size = ws.entities.size();
}

template <typename EvalT>
void assemble_volumetric(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d) {
  auto f = p->get_volumetric();
  for (int block_idx = 0; block_idx < f.size(); ++block_idx) {
    f[block_idx]->template preEvaluate<EvalT>(*(i.get()));
    auto num_ws = d->get_num_worksets(block_idx);
    for (int ws_idx = 0; ws_idx < num_ws; ++ws_idx) {
      load_disc_info(ws, d, block_idx, ws_idx);
      f[block_idx]->template evaluateFields<EvalT>(ws);
    }
    f[block_idx]->template postEvaluate<EvalT>(*(i.get()));
  }
}

template <typename EvalT>
void assemble_neumann(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d) {
  auto f = p->get_neumann();
  f->template preEvaluate<EvalT>(*(i.get()));
  f->template evaluateFields<EvalT>(ws);
  f->template postEvaluate<EvalT>(*(i.get()));
  (void)d;
}

template <typename EvalT>
void assemble_dirichlet(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d) {
  auto f = p->get_dirichlet();
  f->template preEvaluate<EvalT>(*(i.get()));
  f->template evaluateFields<EvalT>(ws);
  f->template postEvaluate<EvalT>(*(i.get()));
  (void)d;
}

void compute_primal_residual(RCP<Physics> p, RCP<SolutionInfo> i,
    RCP<Discretization> d, const double t_current, const double t_previous) {
  typedef goal::Traits::Residual Residual;
  auto t0 = time();
  i->owned->R->putScalar(0.0);
  i->ghost->R->putScalar(0.0);
  Workset ws;
  ws.t_current = t_current;
  ws.t_previous = t_previous;
  assemble_volumetric<Residual>(ws, p, i, d);
  assemble_neumann<Residual>(ws, p, i, d);
  i->gather_R();
  assemble_dirichlet<Residual>(ws, p, i, d);
  auto t1 = time();
  print(" > residual computed in %f seconds", t1 - t0);
}

void compute_primal_jacobian(RCP<Physics> p, RCP<SolutionInfo> i,
    RCP<Discretization> d, const double t_current, const double t_previous) {
  typedef goal::Traits::Jacobian Jacobian;
  auto t0 = time();
  i->owned->R->putScalar(0.0);
  i->ghost->R->putScalar(0.0);
  i->owned->dRdu->resumeFill();
  i->ghost->dRdu->resumeFill();
  i->owned->dRdu->setAllToScalar(0.0);
  i->ghost->dRdu->setAllToScalar(0.0);
  Workset ws;
  ws.t_current = t_current;
  ws.t_previous = t_previous;
  assemble_volumetric<Jacobian>(ws, p, i, d);
  assemble_neumann<Jacobian>(ws, p, i, d);
  i->ghost->dRdu->fillComplete();
  i->gather_R();
  i->gather_dRdu();
  assemble_dirichlet<Jacobian>(ws, p, i, d);
  i->owned->dRdu->fillComplete();
  auto t1 = time();
  print(" > jacobian computed in %f seconds", t1 - t0);
}

void compute_dual_jacobian(RCP<Physics> p, RCP<SolutionInfo> i,
    RCP<Discretization> d, const double t_current, const double t_previous) {
  typedef goal::Traits::Jacobian J;
  auto t0 = time();
  i->owned->dJdu->putScalar(0.0);
  i->ghost->dJdu->putScalar(0.0);
  i->owned->dRdu->resumeFill();
  i->ghost->dRdu->resumeFill();
  i->owned->dRdu->setAllToScalar(0.0);
  i->ghost->dRdu->setAllToScalar(0.0);
  Workset ws;
  ws.t_current = t_current;
  ws.t_previous = t_previous;
  assemble_volumetric<J>(ws, p, i, d);
  assemble_neumann<J>(ws, p, i, d);
  i->gather_R();
  i->ghost->dRdu->fillComplete();
  i->gather_dRdu();
  i->gather_dJdu();
  assemble_dirichlet<J>(ws, p, i, d);
  i->owned->dRdu()->fillComplete();
  auto t1 = time();
  print(" > dual jacobian computed in %f seconds", t1 - t0);
}

void compute_error_residual(RCP<Physics> p, RCP<SolutionInfo> i,
    RCP<Discretization> d, const double t_current, const double t_previous) {
  typedef goal::Traits::Residual R;
  auto t0 = time();
  i->owned->R->putScalar(0.0);
  i->ghost->R->putScalar(0.0);
  Workset ws;
  ws.t_current = t_current;
  ws.t_previous = t_previous;
  assemble_volumetric<R>(ws, p, i, d);
  i->gather_R();
  auto t1 = time();
  print(" > error residual computed in %f seconds", t1 - t0);
}

template void assemble_volumetric<goal::Traits::Residual>(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);
template void assemble_volumetric<goal::Traits::Jacobian>(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);
template void assemble_neumann<goal::Traits::Residual>(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);
template void assemble_neumann<goal::Traits::Jacobian>(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);
template void assemble_dirichlet<goal::Traits::Residual>(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);
template void assemble_dirichlet<goal::Traits::Jacobian>(
    Workset& ws, RCP<Physics> p, RCP<SolutionInfo> i, RCP<Discretization> d);

}  /* namespace goal */
