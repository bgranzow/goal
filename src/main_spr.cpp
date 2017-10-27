#include <Teuchos_YamlParameterListHelpers.hpp>

#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <ma.h>
#include <parma.h>
#include <PCU.h>
#include <spr.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_functional.hpp"
#include "goal_mechanics.hpp"
#include "goal_output.hpp"
#include "goal_primal.hpp"
#include "goal_regression.hpp"
#include "goal_states.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;
using GroupCode = Parma_GroupCode;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<int>("num steps", 0);
  p.set<double>("step size", 0.0);
  p.set<double>("initial time", 0.0);
  p.sublist("discretization");
  p.sublist("dirichlet bcs");
  p.sublist("mechanics");
  p.sublist("functional");
  p.sublist("primal linear algebra");
  p.sublist("adaptation");
  p.sublist("output");
  p.sublist("regression");
  return p;
}

static ParameterList get_valid_adapt_params() {
  ParameterList p;
  p.set<std::string>("field name", "");
  p.set<bool>("is ip field", true);
  p.set<int>("adapt cycles", 1);
  p.set<int>("adapt iters", 1);
  p.set<bool>("fix shape", true);
  p.set<bool>("shoud coarsen", true);
  p.set<double>("good quality", 0.3);
  p.set<int>("target elems", 1);
  p.set<int>("min density", 10000);
  return p;
}

class Solver {
  public:
    Solver(const char* in);
    ~Solver();
    void solve();
  private:
    void adapt(const int step, const int cycle);
    RCP<ParameterList> params;
    Disc* disc;
    Mechanics* mech;
    Primal* primal;
    Functional* functional;
    Output* output;
    int num_cycles;
};

Solver::Solver(const char* in) {
  print("reading input: %s", in);
  params = rcp(new ParameterList);
  Teuchos::updateParametersFromYamlFile(in, params.ptr());
  params->validateParameters(get_valid_params(), 0);
  auto disc_params = params->sublist("discretization");
  auto mech_params = params->sublist("mechanics");
  auto out_params = params->sublist("output");
  auto adapt_params = params->sublist("adaptation");
  disc = create_disc(disc_params);
  mech = create_mechanics(mech_params, disc);
  primal = create_primal(*params, mech);
  functional = create_functional(*params, primal);
  output = create_output(out_params, disc);
  num_cycles = adapt_params.get<int>("adapt cycles");
}

Solver::~Solver() {
  destroy_output(output);
  destroy_functional(functional);
  destroy_primal(primal);
  destroy_mechanics(mech);
  destroy_disc(disc);
}

static double get_avg_density(apf::Mesh* m) {
  double nelem = m->count(m->getDimension());
  nelem = PCU_Add_Double(nelem);
  return nelem / PCU_Comm_Peers();
}

static int get_shrink_factor(apf::Mesh* m, const double min) {
  int factor = 1;
  auto rho = get_avg_density(m);
  while (rho < min) {
    if (factor >= PCU_Comm_Peers()) break;
    factor *= 2;
    rho *= 2;
  }
  GOAL_DEBUG_ASSERT(PCU_Comm_Peers() % factor == 0);
  return factor;
}

static void warn_about_shrink(const int factor) {
  int nprocs = PCU_Comm_Peers() / factor;
  print(" > shrinking to %d procs to avoid turkey gobbling", nprocs);
}

static void run_shrunken(apf::Mesh2* m, const int f, GroupCode& c) {
  if (f == 1) c.run(0);
  else {
    warn_about_shrink(f);
    Parma_ShrinkPartition(m, f, c);
  }
}

static void configure_ma(ma::Input* in, ParameterList& p) {
  in->maximumIterations = p.get<int>("adapt iters", 1);
  in->shouldCoarsen = p.get<bool>("should coarsen", true);
  in->shouldFixShape = p.get<bool>("fix shape", true);
  in->goodQuality = p.get<double>("good quality", 0.2);
}

struct SPRCallback : public GroupCode {
  public:
    SPRCallback(apf::Mesh2* m, const int t, ParameterList const& p) {
      mesh = m;
      target = t;
      params = p;
    };
    void run(int) {
      apf::Field* ip_field = 0;
      auto fname = params.get<std::string>("field name");
      auto is_ip = params.get<bool>("is ip field");
      auto field = mesh->findField(fname.c_str());
      GOAL_DEBUG_ASSERT(field);
      if (is_ip) ip_field = field;
      else ip_field = spr::getGradIPField(field, "spr_grad", 1);
      auto size = spr::getTargetSPRSizeField(ip_field, target);
      if (! is_ip) apf::destroyField(ip_field);
      auto in = ma::configure(mesh, size);
      configure_ma(in, params);
      ma::adapt(in);
      apf::destroyField(size);
    }
  private:
    apf::Mesh2* mesh;
    int target; 
    ParameterList params;
};

void Solver::adapt(const int step, const int cycle) {
  if (cycle == num_cycles) return;
  print("*** adaptation");
  print("*** at step: (%d)", step);
  print("*** and cycle: (%d)", cycle);
  auto mesh = disc->get_apf_mesh();
  auto adapt_params = params->sublist("adaptation");
  auto target = adapt_params.get<int>("target elems");
  auto min_density = adapt_params.get<int>("min density", 10000);
  auto factor = get_shrink_factor(mesh, min_density);
  adapt_params.validateParameters(get_valid_adapt_params(), 0);
  SPRCallback spr(mesh, target, adapt_params);
  run_shrunken(mesh, factor, spr);
}

void Solver::solve() {
  int num_steps = params->get<int>("num steps");
  double dt = params->get<double>("step size");
  double t_old = params->get<double>("initial time");
  double t_now = t_old + dt;
  for (int step = 1; step <= num_steps; ++step) {
    print("****** load step: %d", step);
    print("****** from time: %f", t_old);
    print("****** to time:   %f", t_now);
    for (int cycle = 0; cycle <= num_cycles; ++cycle) {
      disc->build_data();
      primal->build_data();
      primal->solve(t_now, t_old);
      functional->compute(t_now, t_old);
      functional->print_value();
      output->write(t_now, cycle);
      primal->destroy_data();
      disc->destroy_data();
      adapt(step, cycle);
    }
    mech->get_states()->update();
    t_old = t_now;
    t_now += dt;
  }
  check_J_regression(*params, functional);
}

}

int main(int argc, char** argv) {
  goal::initialize();
  GOAL_DEBUG_ASSERT(argc == 2);
  const char* in = argv[1];
  { goal::Solver solver(in);
    solver.solve(); }
  goal::finalize();
}
