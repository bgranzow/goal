#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_assembly.hpp"
#include "goal_avg_disp.hpp"
#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_J2.hpp"
#include "goal_kinematics.hpp"
#include "goal_mechanics.hpp"
#include "goal_mresidual.hpp"
#include "goal_mixed.hpp"
#include "goal_neohookean.hpp"
#include "goal_point_wise.hpp"
#include "goal_presidual.hpp"
#include "goal_scalar_types.hpp"
#include "goal_stabilization.hpp"
#include "goal_states.hpp"

namespace goal {

using Teuchos::rcp;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("model", "");
  p.sublist("materials");
  return p;
}

Mechanics::Mechanics(ParameterList const& p, Disc* d) {
  p.validateParameters(get_valid_params(), 0);
  params = p;
  disc = d;
  displacement = 0;
  pressure = 0;
  states = 0;
  model = params.get<std::string>("model");
  make_displacement();
  make_pressure();
  make_states();
}

Mechanics::~Mechanics() {
  destroy_states(states);
  apf::destroyField(pressure);
  apf::destroyField(displacement);
}

void Mechanics::make_displacement() {
  auto m = disc->get_apf_mesh();
  if (disc->is_base()) {
    displacement = apf::createFieldOn(m, "u", apf::VECTOR);
    apf::zeroField(displacement);
  }
  else {
    displacement = m->findField("u");
    GOAL_DEBUG_ASSERT(displacement);
  }
}

void Mechanics::make_pressure() {
  auto m = disc->get_apf_mesh();
  if (disc->is_base()) {
    pressure = apf::createFieldOn(m, "p", apf::SCALAR);
    apf::zeroField(pressure);
  }
  else {
    pressure = m->findField("p");
    GOAL_DEBUG_ASSERT(pressure);
  }
}

void Mechanics::make_states() {
  states = create_states(disc);
  if (! disc->is_base()) return;
  states->add("sigma", apf::MATRIX);
  if (model == "J2") {
    states->add("eqps", apf::SCALAR, true);
    states->add("Fp", apf::MATRIX, true, true);
  }
}

template <typename T>
void Mechanics::build_resid(Evaluators& E, bool save) {
  ParameterList mat = params.sublist("materials");
  auto u = find_evaluator("u", E);
  auto p = find_evaluator("p", E);
  auto uw = find_evaluator("uw", E);
  auto pw = find_evaluator("pw", E);
  auto kin = rcp(new Kinematics<T>(u));
  RCP<Model<T>> cm;
  if (model == "neohookean")
    cm = rcp(new Neohookean<T>(kin, states, save, mat));
  else if (model == "J2")
    cm = rcp(new J2<T>(kin, states, save, mat));
  else
    fail("unknown model: %s", model.c_str());
  auto mixed = rcp(new Mixed<T>(p, cm, states, save));
  auto mresidual = rcp(new MResidual<T>(u, uw, cm));
  auto presidual = rcp(new PResidual<T>(p, pw, kin, mat));
  auto stab = rcp(new Stabilization<T>(p, pw, cm, kin, mat));
  E.push_back(kin);
  E.push_back(cm);
  E.push_back(mixed);
  E.push_back(mresidual);
  E.push_back(presidual);
  E.push_back(stab);
}

template <typename T>
void Mechanics::build_functional(ParameterList const& params, Evaluators& E) {
  auto type = params.get<std::string>("type");
  auto u = find_evaluator("u", E);
  RCP<QoI<T>> J;
  if (type == "point wise")
    J = rcp(new PointWise<T>(params));
  else if (type == "avg disp")
    J = rcp(new AvgDisp<T>(u));
  else
    fail("unknown functional type: %s", type.c_str());
  E.push_back(J);
}

Mechanics* create_mechanics(ParameterList const& p, Disc* d) {
  return new Mechanics(p, d);
}

void destroy_mechanics(Mechanics* m) {
  delete m;
}

template void Mechanics::build_resid<ST>(Evaluators&, bool);
template void Mechanics::build_resid<FADT>(Evaluators&, bool);
template void Mechanics::build_functional<ST>(ParameterList const&, Evaluators&);
template void Mechanics::build_functional<FADT>(ParameterList const&, Evaluators&);

}
