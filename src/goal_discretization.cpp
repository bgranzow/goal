#include <apf.h>
#include <apfAlbany.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <Teuchos_ParameterList.hpp>

#include "goal_control.hpp"
#include "goal_discretization.hpp"

namespace goal {

using Teuchos::rcp;

static RCP<ParameterList> get_valid_params() {
  auto n = "valid discretization params";
  auto p = rcp(new ParameterList(n));
  p->set<std::string>("geom file", "");
  p->set<std::string>("mesh file", "");
  p->set<std::string>("assoc file", "");
  p->set<bool>("reorder mesh", false);
  p->set<apf::Mesh2*>("mesh", 0);
  p->set<int>("workset size", 0);
  return p;
}

static void validate_params(RCP<const ParameterList> p) {
  assert(p->isType<std::string>("assoc file"));
  assert(p->isType<bool>("reorder mesh"));
  assert(p->isType<int>("workset size"));
  p->validateParameters(*get_valid_params(), 0);
}

static void load_mesh_from_file(apf::Mesh2** mesh, RCP<const ParameterList> p) {
  gmi_register_mesh();
  auto geom_file = p->get<std::string>("geom file");
  auto mesh_file = p->get<std::string>("mesh file");
  auto g = geom_file.c_str();
  auto m = mesh_file.c_str();
  *mesh = apf::loadMdsMesh(g, m);
}

static bool set_mesh(apf::Mesh2** mesh, RCP<const ParameterList> p) {
  bool owns = true;
  if (p->isType<std::string>("geom file") &&
      p->isType<std::string>("mesh file"))
    load_mesh_from_file(mesh, p);
  else if (p->isType<apf::Mesh2*>("mesh")) {
    *mesh = p->get<apf::Mesh2*>("mesh");
    owns = false;
  } else
    fail("unable to set apf mesh");
  return owns;
}

static void reorder_mesh(apf::Mesh2* mesh, RCP<const ParameterList> p) {
  if (p->get<bool>("reorder mesh")) apf::reorderMdsMesh(mesh);
  mesh->verify();
}

static apf::StkModels* read_sets(apf::Mesh* m, RCP<const ParameterList> p) {
  auto sets = new apf::StkModels;
  auto fn = p->get<std::string>("assoc file");
  auto filename = fn.c_str();
  print("reading association file: %s", filename);
  static std::string const setNames[3] = {
      "node set", "side set", "element block"};
  auto d = m->getDimension();
  int dims[3] = {0, d - 1, d};
  std::ifstream f(filename);
  if (!f.good()) fail("cannot open file: %s", filename);
  std::string sline;
  int lc = 0;
  while (std::getline(f, sline)) {
    if (!sline.length()) break;
    ++lc;
    int sdi = -1;
    for (int di = 0; di < 3; ++di)
      if (sline.compare(0, setNames[di].length(), setNames[di]) == 0) sdi = di;
    if (sdi == -1)
      fail("invalid association line # %d:\n\t%s", lc, sline.c_str());
    int sd = dims[sdi];
    std::stringstream strs(sline.substr(setNames[sdi].length()));
    auto set = new apf::StkModel();
    strs >> set->stkName;
    int nents;
    strs >> nents;
    if (!strs) fail("invalid association line # %d:\n\t%s", lc, sline.c_str());
    for (int ei = 0; ei < nents; ++ei) {
      std::string eline;
      std::getline(f, eline);
      if (!f || !eline.length())
        fail("invalid association after line # %d", lc);
      ++lc;
      std::stringstream strs2(eline);
      int mdim, mtag;
      strs2 >> mdim >> mtag;
      if (!strs2) fail("bad associations line # %d:\n\t%s", lc, eline.c_str());
      set->ents.push_back(m->findModelEntity(mdim, mtag));
      if (!set->ents.back())
        fail("no model entity with dim: %d and tag: %d", mdim, mtag);
    }
    sets->models[sd].push_back(set);
  }
  sets->computeInverse();
  return sets;
}

Discretization::Discretization(RCP<const ParameterList> p) {
  params = p;
  validate_params(p);
  owns_mesh = set_mesh(&mesh, params);
  reorder_mesh(mesh, params);
  sets = read_sets(mesh, params);
  num_dims = mesh->getDimension();
  ws_size = params->get<int>("workset size");
  print(" num element blocks: %d", get_num_elem_blocks());
  print(" num side sets:      %d", get_num_side_sets());
  print(" num node sets:      %d", get_num_node_sets());
  update();
}

Discretization::~Discretization() {
  if (mesh && owns_mesh) {
    mesh->destroyNative();
    apf::destroyMesh(mesh);
  }
}

void Discretization::update() {
  double t0 = time();
  compute_elem_blocks();
  compute_side_sets();
  double t1 = time();
  print("discretization updated in %f seconds", t1 - t0);
}

int Discretization::get_num_elem_blocks() const {
  return sets->models[num_dims].size();
}

int Discretization::get_num_side_sets() const {
  return sets->models[num_dims - 1].size();
}

int Discretization::get_num_node_sets() const { return sets->models[0].size(); }

std::string Discretization::get_elem_block_name(const int idx) const {
  assert(idx < get_num_node_sets());
  return sets->models[num_dims][idx]->stkName;
}

std::string Discretization::get_side_set_name(const int idx) const {
  assert(idx < get_num_side_sets());
  return sets->models[num_dims - 1][idx]->stkName;
}

std::string Discretization::get_node_set_name(const int idx) const {
  assert(idx < get_num_node_sets());
  return sets->models[0][idx]->stkName;
}

int Discretization::get_elem_type(const int block_idx) {
  auto name = get_elem_block_name(block_idx);
  auto e = elem_blocks[name][0][0];
  return mesh->getType(e);
}

int Discretization::get_num_worksets(const int block_idx) {
  auto name = get_elem_block_name(block_idx);
  return elem_blocks[name].size();
}

std::vector<apf::MeshEntity*> const& Discretization::get_elems(
    std::string const& elem_block, const int ws_idx) {
  assert(elem_blocks.count(elem_block));
  return elem_blocks[elem_block][ws_idx];
}

std::vector<apf::MeshEntity*> const& Discretization::get_sides(
    std::string const& side_set) {
  assert(side_sets.count(side_set));
  return side_sets[side_set];
}

void Discretization::compute_elem_blocks() {
  int neb = get_num_elem_blocks();
  for (int i = 0; i < neb; ++i) elem_blocks[get_elem_block_name(i)].resize(0);
  apf::MeshEntity* elem;
  auto it = mesh->begin(num_dims);
  std::map<std::string, std::vector<apf::MeshEntity*> > map;
  for (int i = 0; i < neb; ++i) map[get_elem_block_name(i)].resize(0);
  while ((elem = mesh->iterate(it))) {
    auto mr = mesh->toModel(elem);
    auto stkm = sets->invMaps[num_dims][mr];
    auto name = stkm->stkName;
    if (map[name].size() >= std::size_t(ws_size)) {
      elem_blocks[name].push_back(map[name]);
      map[name].clear();
    }
    map[name].push_back(elem);
  }
  mesh->end(it);
  for (int i = 0; i < neb; ++i) {
    auto name = get_elem_block_name(i);
    elem_blocks[name].push_back(map[name]);
  }
}

void Discretization::compute_side_sets() {
  int nss = get_num_side_sets();
  for (int i = 0; i < nss; ++i) side_sets[get_side_set_name(i)].resize(0);
  apf::MeshEntity* side;
  auto it = mesh->begin(num_dims - 1);
  while ((side = mesh->iterate(it))) {
    auto me = mesh->toModel(side);
    if (!sets->invMaps[num_dims - 1].count(me)) continue;
    auto ss = sets->invMaps[num_dims - 1][me];
    auto ssn = ss->stkName;
    apf::Up adj_elems;
    mesh->getUp(side, adj_elems);
    assert(adj_elems.n == 1);
    side_sets[ssn].push_back(side);
  }
  mesh->end(it);
}

}  // namespace goal
