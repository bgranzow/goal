#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_assembly.hpp"
#include "goal_avg_disp.hpp"
#include "goal_avg_disp_subdomain.hpp"
#include "goal_avg_vm.hpp"
#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_J2.hpp"
#include "goal_kinematics.hpp"
#include "goal_ks_vm.hpp"
#include "goal_mechanics.hpp"
#include "goal_mresidual.hpp"
#include "goal_mixed.hpp"
#include "goal_neohookean.hpp"
#include "goal_point_wise.hpp"
#include "goal_presidual.hpp"
#include "goal_scalar_types.hpp"
#include "goal_stabilization.hpp"
#include "goal_states.hpp"
#include "goal_temp.hpp"

namespace goal {

using Teuchos::rcp;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("model", "");
  p.sublist("temperature");
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
  have_temp = false;
  if (params.isSublist("temperature")) {
    temp_params = params.sublist("temperature");
    have_temp = true;
  }
}

Mechanics::~Mechanics() {
  destroy_states(states);
  apf::destroyField(pressure);
  apf::destroyField(displacement);
}

void Mechanics::make_displacement() {
  auto m = disc->get_apf_mesh();
  auto f = m->findField("u");
  if (f) displacement = f;
  else displacement = apf::createFieldOn(m, "u", apf::VECTOR);
  if (!f) apf::zeroField(displacement);
}

void Mechanics::make_pressure() {
  auto m = disc->get_apf_mesh();
  auto f = m->findField("p");
  if (f) pressure = f;
  else pressure = apf::createFieldOn(m, "p", apf::SCALAR);
  if (!f) apf::zeroField(pressure);
}

void Mechanics::make_states() {
  states = create_states(disc);
  states->add("sigma", apf::MATRIX);
  if (model == "J2") {
    states->add("eqps", apf::SCALAR, true);
    states->add("Fp", apf::MATRIX, true, true);
  }
  states->update();
}

template <typename T>
void Mechanics::build_resid(Evaluators& E, bool save) {

  ParameterList mat = params.sublist("materials");

  auto u = find_evaluator("u", E);
  auto p = find_evaluator("p", E);
  auto uw = find_evaluator("uw", E);
  auto pw = find_evaluator("pw", E);

  auto kin = rcp(new Kinematics<T>(u));
  E.push_back(kin);

  RCP<Model<T>> cm;
  if (model == "neohookean")
    cm = rcp(new Neohookean<T>(kin, states, save, mat));
  else if (model == "J2")
    cm = rcp(new J2<T>(kin, states, save, mat));
  else
    fail("unknown model: %s", model.c_str());
  E.push_back(cm);

  auto mixed = rcp(new Mixed<T>(p, cm, states, save));
  E.push_back(mixed);

  if (have_temp) {
    auto temp = rcp(new Temp<T>(temp_params, mat, cm, kin, states, save));
    E.push_back(temp);
  }

  auto mresidual = rcp(new MResidual<T>(u, uw, cm));
  E.push_back(mresidual);

  auto presidual = rcp(new PResidual<T>(p, pw, kin, mat));
  E.push_back(presidual);

  auto stab = rcp(new Stabilization<T>(p, pw, kin, mat));
  E.push_back(stab);
}

template <typename T>
void Mechanics::build_functional(ParameterList const& params, Evaluators& E) {
  auto type = params.get<std::string>("type");
  auto u = find_evaluator("u", E);
  auto model = find_evaluator("model", E);
  RCP<QoI<T>> J;
  if (type == "point wise")
    J = rcp(new PointWise<T>(params));
  else if (type == "avg disp")
    J = rcp(new AvgDisp<T>(u));
  else if (type == "avg disp subdomain")
    J = rcp(new AvgDispSubdomain<T>(params, u));
  else if (type == "avg vm")
    J = rcp(new AvgVM<T>(params, model));
  else if (type == "max vm")
    J = rcp(new KSVM<T>(params, model));
  else
    fail("unknown functional type: %s", type.c_str());
  E.push_back(J);
}

void Mechanics::build_error(Evaluators& E) {

  ParameterList mat = params.sublist("materials");

  auto u = find_evaluator("u", E);
  auto p = find_evaluator("p", E);
  auto uw = find_evaluator("uw", E);
  auto pw = find_evaluator("pw", E);
  auto pwc = find_evaluator("pwc", E);

  auto kin = rcp(new Kinematics<ST>(u));
  E.push_back(kin);

  RCP<Model<ST>> cm;
  if (model == "neohookean")
    cm = rcp(new Neohookean<ST>(kin, states, false, mat));
  else if (model == "J2")
    cm = rcp(new J2<ST>(kin, states, false, mat));
  else
    fail("unknown model: %s", model.c_str());
  E.push_back(cm);

  auto mixed = rcp(new Mixed<ST>(p, cm, states, false));
  E.push_back(mixed);

  if (have_temp) {
    auto temp = rcp(new Temp<ST>(temp_params, mat, cm, kin, states, false));
    E.push_back(temp);
  }

  auto mresidual = rcp(new MResidual<ST>(u, uw, cm));
  E.push_back(mresidual);

  auto presidual = rcp(new PResidual<ST>(p, pw, kin, mat));
  E.push_back(presidual);

  auto stab = rcp(new Stabilization<ST>(p, pwc, kin, mat));
  E.push_back(stab);
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
