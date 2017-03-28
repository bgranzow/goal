/** [Step2 all] */
#include <apf.h>
#include <apfBox.h>
#include <apfMesh2.h>
#include <apfAlbany.h>
#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_output.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Step2 {

/** [Step2 using] */
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;

/** [Step2 using] */

/** [Step2 initialize] */
static apf::Mesh2* init() {
  goal::initialize();
  goal::print("Goal: Step 2");
  auto mesh = apf::makeMdsBox(20, 20, 0, 1.0, 1.0, 0.0, true);
  mesh->verify();
  return mesh;
}

/** [Step2 initialize] */

/** [Step2 create_assoc] */
static apf::StkModels* create_associations(apf::Mesh2* m) {
  int dim = m->getDimension();
  auto sets = new apf::StkModels;
  auto square = new apf::StkModel;
  square->stkName = "square";
  square->ents.push_back(m->findModelEntity(2, 0));
  sets->models[dim].push_back(square);
  auto edges = new apf::StkModel;
  edges->stkName = "edges";
  for (int i = 0; i < 4; ++i)
    edges->ents.push_back(m->findModelEntity(1, i));
  sets->models[0].push_back(edges);
  sets->computeInverse();
  return sets;
}

/** [Step2 create_assoc] */

/** [Step2 create_disc] */
static RCP<goal::Discretization> create_disc(
    apf::Mesh2* m, apf::StkModels* s) {
  auto p = rcp(new ParameterList);
  p->set<apf::Mesh2*>("mesh", m);
  p->set<apf::StkModels*>("associations", s);
  p->set<bool>("reorder mesh", true);
  p->set<int>("workset size", 1000);
  return rcp(new goal::Discretization(p));
}

/** [Step2 create_disc] */

/** [Step2 output] */
static void output(RCP<goal::Discretization> d) {
  auto p = rcp(new ParameterList);
  p->set<std::string>("out file", "out_step2");
  auto output = rcp(new goal::Output(p, d));
  output->write(0.0);
}

/** [Step2 output] */

/** [Step2 clean_up] */
static void clean_up(
    apf::Mesh2* m, apf::StkModels* s, RCP<goal::Discretization> d) {
  d = Teuchos::null;
  delete s;
  apf::destroyMesh(m);
}

/** [Step2 clean_up] */

}  // namespace Step2

/** [Step2 main] */
int main() {
  auto mesh = Step2::init();
  auto sets = Step2::create_associations(mesh);
  auto disc = Step2::create_disc(mesh, sets);
  Step2::output(disc);
  Step2::clean_up(mesh, sets, disc);
}

/** [Step2 main] */
/** [Step2 all] */
