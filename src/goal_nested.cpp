#include "goal_control.hpp"
#include "goal_nested.hpp"

#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
#include <PCU.h>

namespace goal {

Nested::Nested(Disc* d, int m) {
  double t0 = time();
  mode = m;
  ratio = 0.0;
  am_base = false;
  sets = d->get_model_sets();
  base_mesh = d->get_apf_mesh();
  number_elems();
  copy_mesh();
  tag_old_verts();
  refine_mesh();
  store_old_verts();
  initialize();
  double t1 = time();
  print(" > nested mesh built in %f seconds", t1 - t0);
}

Nested::~Nested() {
  auto e1 = base_mesh->findField("elems");
  auto e2 = mesh->findField("elems");
  apf::destroyField(e1);
  apf::destroyField(e2);
  apf::removeTagFromDimension(mesh, old_vtx_tag, 0);
  apf::removeTagFromDimension(mesh, new_vtx_tag, 0);
  mesh->destroyTag(old_vtx_tag);
  mesh->destroyTag(new_vtx_tag);
}

void Nested::number_elems() {
  int i = 0;
  int dim = base_mesh->getDimension();
  base_elems.resize(base_mesh->count(dim));
  auto vt = apf::SCALAR;
  auto s = apf::getIPFitShape(dim, 1);
  auto f = apf::createField(base_mesh, "elems", vt, s);
  apf::MeshEntity* elem;
  auto it = base_mesh->begin(dim);
  while ((elem = base_mesh->iterate(it))) {
    base_elems[i] = elem;
    apf::setScalar(f, elem, 0, i++);
  }
  base_mesh->end(it);
}

void Nested::copy_mesh() {
  auto model = base_mesh->getModel();
  mesh = apf::createMdsMesh(model, base_mesh);
  apf::disownMdsModel(mesh);
}

void Nested::tag_old_verts() {
  int i = -1;
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = mesh->begin(0);
  old_vtx_tag = mesh->createIntTag("ovt", 1);
  new_vtx_tag = mesh->createIntTag("nvt", 2);
  while ((vtx = mesh->iterate(it)))
    mesh->setIntTag(vtx, old_vtx_tag, &(++i));
  mesh->end(it);
  old_vertices.resize(mesh->count(0));
}

class Transfer : public ma::SolutionTransfer {
  public:
    Transfer(apf::Mesh* m, apf::MeshTag* o, apf::MeshTag* n);
    bool hasNodesOn(int dim);
    void onVertex(apf::MeshElement* p, ma::Vector&, ma::Entity* vtx);
    void onRefine(ma::Entity* p, ma::EntityArray& c);
    int num_elems;
    double ratio;
  private:
    apf::Mesh* mesh;
    apf::MeshTag* old_vtx_tag;
    apf::MeshTag* new_vtx_tag;
    int num_dims;
};

Transfer::Transfer(apf::Mesh* m, apf::MeshTag* o, apf::MeshTag* n) {
  mesh = m;
  old_vtx_tag = o;
  new_vtx_tag = n;
  num_dims = mesh->getDimension();
  num_elems = 0;
  ratio = 0.0;
}

bool Transfer::hasNodesOn(int dim) {
  if (dim == 0) return true;
  if (dim == num_dims) return true;
  return false;
}

void Transfer::onVertex(apf::MeshElement* p, ma::Vector&, ma::Entity* vtx) {
  int tags[2];
  apf::MeshEntity* verts[2];
  auto edge = apf::getMeshEntity(p);
  mesh->getDownward(edge, 0, verts);
  mesh->getIntTag(verts[0], old_vtx_tag, &(tags[0]));
  mesh->getIntTag(verts[1], old_vtx_tag, &(tags[1]));
  mesh->setIntTag(vtx, new_vtx_tag, &(tags[0]));
}

static double get_size(apf::Mesh* m, apf::MeshEntity* e) {
  double h = 0.0;
  apf::Downward edges;
  int ne = m->getDownward(e, 1, edges);
  for (int i = 0; i < ne; ++i)
    h += apf::measure(m, edges[i]) * apf::measure(m, edges[i]);
  return std::sqrt(h/ne);
}

void Transfer::onRefine(ma::Entity* p, ma::EntityArray& c) {
  auto h_old = get_size(mesh, p);
  for (size_t i = 0; i < c.getSize(); ++i) {
    double h_new = get_size(mesh, c[i]);
    ratio += h_new / h_old;
    num_elems++;
  }
}

static double get_ratio(Transfer* t) {
  auto ratio = t->ratio;
  auto num_elems = t->num_elems;
  PCU_Add_Doubles(&ratio, 1);
  PCU_Add_Ints(&num_elems, 1);
  return ratio / num_elems;
}

void Nested::refine_uniform() {
  goal::print(" > nested: uniform");
  ma::AutoSolutionTransfer transfers(mesh);
  auto mytransfer = new Transfer(mesh, old_vtx_tag, new_vtx_tag);
  transfers.add(mytransfer);
  auto in = ma::configureUniformRefine(mesh, 1, &transfers);
  in->shouldFixShape = false;
  in->shouldSnap = false;
  ma::adapt(in);
  ratio = get_ratio(mytransfer);
}

struct Long : public ma::IdentitySizeField {
  Long(apf::Mesh2* m);
  bool shouldSplit(apf::MeshEntity* edge);
  void indicate_edges();
  void gather_edges();
  void mark_edges();
  int dim;
  apf::Mesh2* mesh;
  apf::Field* marks;
  std::set<apf::MeshEntity*> stored;
};

Long::Long(apf::Mesh2* m) : ma::IdentitySizeField(m) {
  mesh = m;
  dim = mesh->getDimension();
  mark_edges();
}

bool Long::shouldSplit(apf::MeshEntity* edge) {
  if (stored.count(edge)) return true;
  else return false;
}

static apf::MeshEntity* get_max(apf::Mesh* m, apf::MeshEntity* e) {
  apf::Downward edges;
  apf::MeshEntity* max = 0;
  int ne = m->getDownward(e, 1, edges);
  for (int i = 0; i < ne-1; ++i) {
    double hi = apf::measure(m, edges[i]);
    double hip = apf::measure(m, edges[i+1]);
    if (hi > hip) max = edges[i];
    else max = edges[i+1];
  }
  return max;
}

void Long::indicate_edges() {
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = mesh->begin(dim);
  while ((elem = mesh->iterate(elems))) {
    auto edge = get_max(mesh, elem);
    apf::setScalar(marks, edge, 0, 1.0);
  }
  mesh->end(elems);
  apf::accumulate(marks);
}

void Long::gather_edges() {
  apf::MeshEntity* edge;
  apf::MeshIterator* edges = mesh->begin(1);
  while ((edge = mesh->iterate(edges)))
    if (apf::getScalar(marks, edge, 0) > 0)
      stored.insert(edge);
  mesh->end(edges);
}

void Long::mark_edges() {
  auto s = apf::getConstant(1);
  marks = apf::createField(mesh, "mark", apf::SCALAR, s);
  apf::zeroField(marks);
  indicate_edges();
  gather_edges();
  apf::destroyField(marks);
}

void Nested::refine_long() {
  goal::print(" > nested: long");
  ma::AutoSolutionTransfer transfers(mesh);
  auto mytransfer = new Transfer(mesh, old_vtx_tag, new_vtx_tag);
  transfers.add(mytransfer);
  auto size = new Long(mesh);
  auto in = ma::configureIdentity(mesh, size, &transfers);
  in->shouldFixShape = false;
  in->shouldSnap = false;
  in->maximumIterations = 1;
  ma::adapt(in);
  delete size;
}

void Nested::refine_mesh() {
  if (mode == UNIFORM) refine_uniform();
  else if (mode == LONG) refine_long();
  else fail("unknown refine mode: %d", mode);
  apf::reorderMdsMesh(mesh);
  mesh->verify();
}

void Nested::store_old_verts() {
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = mesh->begin(0);
  while ((vtx = mesh->iterate(it))) {
    if (! mesh->hasTag(vtx, old_vtx_tag)) continue;
    int tag;
    mesh->getIntTag(vtx, old_vtx_tag, &tag);
    old_vertices[tag] = vtx;
  }
  mesh->end(it);
}

Nested* create_nested(Disc* d, int m) {
  return new Nested(d, m);
}

void destroy_nested(Nested* n) {
  delete n;
}

}
