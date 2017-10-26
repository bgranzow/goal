#include <apf.h>
#include <apfMesh2.h>

#include "goal_assembly.hpp"
#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_integrator.hpp"
#include "goal_sol_info.hpp"

namespace goal {

RCP<Integrator> find_evaluator(std::string const& n, Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    if (E[i]->get_name() == n)
      return E[i];
  return Teuchos::null;
}

void set_time(Evaluators& E, const double t_now, const double t_old) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->set_time(t_now, t_old);
}

void pre_process(SolInfo* s, Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->pre_process(s);
}

void set_elem_sets(const int es_idx, Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->set_elem_set(es_idx);
}

void gather(apf::MeshElement* me, Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->gather(me);
}

void in_elem(apf::MeshElement* me, Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->in_elem(me);
}

void at_point(
    apf::Vector3 const& p, double w, double dv, Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->at_point(p, w, dv);
}

void out_elem(Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->out_elem();
}

void scatter(SolInfo* s, Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->scatter(s);
}

void post_process(SolInfo* s, Evaluators const& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->post_process(s);
}

void assemble(Evaluators const& E, SolInfo* s) {
  apf::Vector3 xi;
  auto disc = s->get_disc();
  auto mesh = disc->get_apf_mesh();
  pre_process(s, E);
  for (int es = 0; es < disc->get_num_elem_sets(); ++es) {
    set_elem_sets(es, E);
    auto esn = disc->get_elem_set_name(es);
    auto elems = disc->get_elems(esn);
    for (size_t elem = 0; elem < elems.size(); ++elem) {
      auto me = apf::createMeshElement(mesh, elems[elem]);
      gather(me, E);
      in_elem(me, E);
      apf::getIntPoint(me, 1, 0, xi);
      double dv = apf::getDV(me, xi);
      double w = apf::getIntWeight(me, 1, 0);
      at_point(xi, w, dv, E);
      out_elem(E);
      scatter(s, E);
      apf::destroyMeshElement(me);
    }
  }
  post_process(s, E);
}

}
