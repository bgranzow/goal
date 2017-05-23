#include <apf.h>
#include <apfAlbany.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <Teuchos_ParameterList.hpp>

#include "goal_control.hpp"
#include "goal_discretization.hpp"

namespace goal {

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("geom file", "");
  p.set<std::string>("mesh file", "");
  p.set<std::string>("assoc file", "");
  p.set<bool>("reorder mesh", "");
  p.set<bool>("make quadratic", "");
  p.set<int>("workset size", 0);
  p.set<apf::Mesh2*>("mesh", 0);
  p.set<apf::StkModels*>("associations", 0);
  return p;
}

static void validate_params(ParameterList const& p) {
  GOAL_ALWAYS_ASSERT_VERBOSE(p.isType<bool>("reorder mesh"),
      "'reorder mesh' is not a discretization parameter");
  GOAL_ALWAYS_ASSERT_VERBOSE(p.isType<int>("workset size"),
      "'workset size' is not a discretization parameter");
  p.validateParameters(get_valid_params(), 0);
}

static void load_mesh_from_file(apf::Mesh2** mesh, ParameterList const& p) {
  gmi_register_mesh();
  auto geom_file = p.get<std::string>("geom file");
  auto mesh_file = p.get<std::string>("mesh file");
  auto g = geom_file.c_str();
  auto m = mesh_file.c_str();
  *mesh = apf::loadMdsMesh(g, m);
}

static bool set_mesh(apf::Mesh2** mesh, ParameterList const& p) {
  bool owns = true;
  if (p.isType<std::string>("geom file") &&
      p.isType<std::string>("mesh file"))
    load_mesh_from_file(mesh, p);
  else if (p.isType<apf::Mesh2*>("mesh")) {
    *mesh = p.get<apf::Mesh2*>("mesh");
    owns = false;
  } else
    fail("unable to set apf mesh");
  return owns;
}

static void make_quadratic(apf::Mesh2* mesh, ParameterList const& p) {
  if (! p.isType<bool>("make quadratic")) return;
  if (! p.get<bool>("make quadratic")) return;
  if (mesh->getShape()->getOrder() == 2) return;
  apf::changeMeshShape(mesh, apf::getSerendipity());
}

static void reorder_mesh(apf::Mesh2* mesh, ParameterList const& p) {
  if (p.get<bool>("reorder mesh"))
    apf::reorderMdsMesh(mesh);
}

static apf::StkModels* read_sets(apf::Mesh* m, ParameterList const& p) {
  auto sets = new apf::StkModels;
  auto fn = p.get<std::string>("assoc file");
  auto filename = fn.c_str();
  print("reading association file: %s", filename);
  static std::string const setNames[3] = {
      "node set", "side set", "element set"};
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

static bool set_associations(
    apf::StkModels** sets, apf::Mesh2* m, ParameterList const& p) {
  bool owns = true;
  if (p.isType<std::string>("assoc file"))
    *sets = read_sets(m, p);
  else if (p.isType<apf::StkModels*>("associations")) {
    *sets = p.get<apf::StkModels*>("associations");
    owns = false;
  } else
    fail("unable to set apf associations");
  return owns;
}

Discretization::Discretization(ParameterList const& p) {
  params = p;
  validate_params(p);
  owns_mesh = set_mesh(&mesh, params);
  make_quadratic(mesh, params);
  reorder_mesh(mesh, params);
  mesh->verify();
  owns_sets = set_associations(&sets, mesh, params);
  num_dims = mesh->getDimension();
  ws_size = params.get<int>("workset size");
  print(" num element sets:   %d", get_num_elem_sets());
  print(" num side sets:      %d", get_num_side_sets());
  print(" num node sets:      %d", get_num_node_sets());
  update();
}

Discretization::~Discretization() {
  if (mesh && owns_mesh) {
    mesh->destroyNative();
    apf::destroyMesh(mesh);
  }
  if (sets && owns_sets)
    delete sets;
}

void Discretization::update() {
  double t0 = time();
  compute_elem_sets();
  compute_side_sets();
  double t1 = time();
  print("discretization updated in %f seconds", t1 - t0);
}

int Discretization::get_num_elem_sets() const {
  return sets->models[num_dims].size();
}

int Discretization::get_num_side_sets() const {
  return sets->models[num_dims - 1].size();
}

int Discretization::get_num_node_sets() const {
  return sets->models[0].size();
}

std::string Discretization::get_elem_set_name(const int i) const {
  GOAL_DEBUG_ASSERT(i < get_num_elem_sets());
  return sets->models[num_dims][i]->stkName;
}

std::string Discretization::get_side_set_name(const int i) const {
  GOAL_DEBUG_ASSERT(i < get_num_side_sets());
  return sets->models[num_dims - 1][i]->stkName;
}

std::string Discretization::get_node_set_name(const int i) const {
  GOAL_DEBUG_ASSERT(i < get_num_node_sets());
  return sets->models[0][i]->stkName;
}

int Discretization::get_elem_set_idx(std::string const& n) const {
  int idx = -1;
  for (int i = 0; i < get_num_elem_sets(); ++i)
    if (n == get_elem_set_name(i))
      idx = i;
  return idx;
}

int Discretization::get_side_set_idx(std::string const& n) const {
  int idx = -1;
  for (int i = 0; i < get_num_side_sets(); ++i)
    if (n == get_side_set_name(i))
      idx = i;
  return idx;
}

int Discretization::get_node_set_idx(std::string const& n) const {
  int idx = -1;
  for (int i = 0; i < get_num_node_sets(); ++i)
    if (n == get_node_set_name(i))
      idx = i;
  return idx;
}

int Discretization::get_elem_type(const int i) {
  auto name = get_elem_set_name(i);
  if (elem_sets[name].size() == 0)
    return -1;
  auto e = elem_sets[name][0][0];
  return mesh->getType(e);
}

int Discretization::get_side_type(const int i) {
  auto name = get_side_set_name(i);
  if (side_sets[name].size() == 0)
    return -1;
  auto s = side_sets[name][0][0];
  return mesh->getType(s);
}

int Discretization::get_num_elem_worksets(const int i) {
  auto name = get_elem_set_name(i);
  return elem_sets[name].size();
}

int Discretization::get_num_side_worksets(const int i) {
  auto name = get_side_set_name(i);
  return side_sets[name].size();
}

std::vector<apf::MeshEntity*> const& Discretization::get_elems(
    std::string const& elem_set, const int ws_idx) {
  if (! elem_sets.count(elem_set))
    fail("element set %s not found", elem_set.c_str());
  return elem_sets[elem_set][ws_idx];
}

std::vector<apf::MeshEntity*> const& Discretization::get_sides(
    std::string const& side_set, const int ws_idx) {
  if (! side_sets.count(side_set))
    fail("side set %s not found", side_set.c_str());
  return side_sets[side_set][ws_idx];
}

void Discretization::compute_elem_sets() {
  int neb = get_num_elem_sets();
  for (int i = 0; i < neb; ++i)
    elem_sets[ get_elem_set_name(i) ].resize(0);
  std::map<std::string, std::vector<apf::MeshEntity*> > map;
  for (int i = 0; i < neb; ++i)
    map[ get_elem_set_name(i) ].resize(0);
  apf::MeshEntity* elem;
  auto it = mesh->begin(num_dims);
  while ((elem = mesh->iterate(it))) {
    auto mr = mesh->toModel(elem);
    auto stkm = sets->invMaps[num_dims][mr];
    auto name = stkm->stkName;
    if (map[name].size() >= std::size_t(ws_size)) {
      elem_sets[name].push_back(map[name]);
      map[name].clear();
    }
    map[name].push_back(elem);
  }
  mesh->end(it);
  for (int i = 0; i < neb; ++i) {
    auto name = get_elem_set_name(i);
    if (map[name].size() > 0)
      elem_sets[name].push_back(map[name]);
  }
}

void Discretization::compute_side_sets() {
  int nss = get_num_side_sets();
  for (int i = 0; i < nss; ++i)
    side_sets[ get_side_set_name(i) ].resize(0);
  std::map<std::string, std::vector<apf::MeshEntity*> > map;
  for (int i = 0; i < nss; ++i)
    map[ get_side_set_name(i) ].resize(0);
  apf::MeshEntity* side;
  auto it = mesh->begin(num_dims - 1);
  while ((side = mesh->iterate(it))) {
    auto me = mesh->toModel(side);
    if (!sets->invMaps[num_dims - 1].count(me)) continue;
    auto stkm = sets->invMaps[num_dims -1][me];
    auto name = stkm->stkName;
    apf::Up adj_elems;
    mesh->getUp(side, adj_elems);
    GOAL_DEBUG_ASSERT_VERBOSE(adj_elems.n == 1,
        "side set defined on non-manifold geometric entity?");
    if (map[name].size() >= std::size_t(ws_size)) {
      side_sets[name].push_back(map[name]);
      map[name].clear();
    }
    map[name].push_back(side);
  }
  mesh->end(it);
  for (int i = 0; i < nss; ++i) {
    auto name = get_side_set_name(i);
    if (map[name].size() > 0)
      side_sets[name].push_back(map[name]);
  }
}

Discretization* create_disc(ParameterList const& p) {
  return new Discretization(p);
}

void destroy_disc(Discretization* d) {
  delete d;
}

} // end namespace goal
