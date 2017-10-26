#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <PCU.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_output.hpp"

namespace goal {

using Teuchos::Array;

static ParameterList get_valid_params() {
  ParameterList p;
  Array<std::string> dummy(0);
  p.set<std::string>("out file", "");
  p.set<int>("interval", 0);
  p.set<bool>("turn off", false);
  return p;
}

static void write_initial_pvd(std::string const& n, int& pos) {
  if (PCU_Comm_Self()) return;
  auto pvd = n + ".pvd";
  std::fstream pvdf;
  pvdf.open(pvd.c_str(), std::ios::out);
  pvdf << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
  pvdf << "  <Collection>" << std::endl;
  pos = pvdf.tellp();
  pvdf << "  </Collection>" << std::endl;
  pvdf << "</VTKFile>" << std::endl;
  pvdf.close();
}

Output::Output(ParameterList const& p, Disc* d) :
    disc(d),
    params(p),
    interval(1),
    turn_off(false),
    pos(0),
    index(0) {
  params.validateParameters(get_valid_params(), 0);
  name = params.get<std::string>("out file");
  if (params.isParameter("turn off"))
    turn_off = params.get<bool>("turn off");
  if (params.isParameter("interval"))
    interval = params.get<int>("interval");
  if (! turn_off) write_initial_pvd(name, pos);
}

Output::~Output() {
}

static void update_pvd(
    std::string const& name,
    std::string const& vtu,
    int& pos,
    const double t) {
  if (PCU_Comm_Self()) return;
  std::string pvd = name + ".pvd";
  std::fstream pvdf;
  pvdf.open(pvd.c_str(), std::ios::out | std::ios::in);
  pvdf.seekp(pos);
  pvdf << "    <DataSet timestep=\"" << t << "\" group=\"\" ";
  pvdf << "part=\"0\" file=\"" << vtu << "/" << vtu;
  pvdf << ".pvtu\"/>" << std::endl;
  pos = pvdf.tellp();
  pvdf << "  </Collection>" << std::endl;
  pvdf << "</VTKFile>" << std::endl;
  pvdf.close();
}

void Output::write_vtk(const double t) {
  auto m = disc->get_apf_mesh();
  std::ostringstream oss;
  oss << name << "_" << index;
  std::string vtu = oss.str();
  update_pvd(name, vtu, pos, t);
  apf::writeVtkFiles(vtu.c_str(), m);
  ++index;
}

void Output::write(const double t, const int iter) {
  if (turn_off) return;
  static int my_out_interval = 0;
  if (my_out_interval++ % interval) return;
  double eps = 1.0e-4;
  write_vtk(t + eps*iter);
}

Output* create_output(ParameterList const& p, Disc* d) {
  return new Output(p, d);
}

void destroy_output(Output* o) {
  delete o;
}

}
