#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_scalar_types.hpp"
#include "goal_states.hpp"

namespace goal {

static ST get_val(ST const& v) { return v; }
static ST get_val(FADT const& v) { return v.val(); }

static void zero(apf::Matrix3x3& v) {
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j)
    v[i][j] = 0.0;
}

template <typename T>
void get_scalar(Disc* d, const char* n, apf::MeshEntity* e, T& v) {
  auto m = d->get_apf_mesh();
  auto f = m->findField(n);
  auto val = apf::getScalar(f, e, 0);
  v = (T)val;
}

template <typename T>
void get_tensor(Disc* d, const char* n, apf::MeshEntity* e, minitensor::Tensor<T>& v) {
  apf::Matrix3x3 val;
  auto m = d->get_apf_mesh();
  auto f = m->findField(n);
  apf::getMatrix(f, e, 0, val);
  for (size_t i = 0; i < v.get_dimension(); ++i)
  for (size_t j = 0; j < v.get_dimension(); ++j)
    v(i,j) = (T)val[i][j];
}

template <typename T>
void set_scalar(Disc* d, const char* n, apf::MeshEntity* e, T const& v) {
  auto m = d->get_apf_mesh();
  auto f = m->findField(n);
  apf::setScalar(f, e, 0, get_val(v));
}

template <typename T>
void set_tensor(Disc* d, const char* n, apf::MeshEntity* e, minitensor::Tensor<T> const& v) {
  apf::Matrix3x3 val;
  zero(val);
  for (size_t i = 0; i < v.get_dimension(); ++i)
  for (size_t j = 0; j < v.get_dimension(); ++j)
    val[i][j] = get_val(v(i,j));
  auto m = d->get_apf_mesh();
  auto f = m->findField(n);
  apf::setMatrix(f, e, 0, val);
}

States::States(Disc* d) {
  disc = d;
  mesh = disc->get_apf_mesh();
  num_dims = disc->get_num_dims();
}

States::~States() {
  for (size_t i = 0; i < states.size(); ++i)
    apf::destroyField(states[i]);
  for (size_t i = 0; i < old_states.size(); ++i)
    apf::destroyField(old_states[i]);
}

void States::add_state(const char* n, int type) {
  auto s = apf::getIPFitShape(num_dims, 1);
  auto f = mesh->findField(n);
  if (!f) f = apf::createField(mesh, n, type, s);
  states.push_back(f);
}

void States::add_old_state(const char* n, int type) {
  auto s = apf::getIPFitShape(num_dims, 1);
  auto on = (std::string)n + "_old";
  auto f = mesh->findField(on.c_str());
  if (!f) f = apf::createField(mesh, on.c_str(), type, s);
  old_states.push_back(f);
}

static apf::Matrix3x3 get_I(int d) {
  static apf::Matrix3x3 I1x1(1,0,0,0,0,0,0,0,0);
  static apf::Matrix3x3 I2x2(1,0,0,0,1,0,0,0,0);
  static apf::Matrix3x3 I3x3(1,0,0,0,1,0,0,0,1);
  switch (d) {
    case 1: return I1x1; break;
    case 2: return I2x2; break;
    case 3: return I3x3; break;
    default: fail("unable to set identity");
  }
}

static void set_I(apf::Field* f) {
  GOAL_ALWAYS_ASSERT(apf::getValueType(f) == apf::MATRIX);
  auto m = apf::getMesh(f);
  auto d = m->getDimension();
  auto I = get_I(d);
  apf::MeshEntity* elem;
  auto elems = m->begin(d);
  while ((elem = m->iterate(elems))) {
    auto type = m->getType(elem);
    auto nodes = apf::getShape(f)->countNodesOn(type);
    for (int n = 0; n < nodes; ++n)
      apf::setMatrix(f, elem, n, I);
  }
  m->end(elems);
}

static void init_state(apf::Field* f, bool I) {
  if (I) set_I(f);
  else apf::zeroField(f);
}

void States::add(const char* n, int t, bool s, bool I) {
  add_state(n, t);
  init_state(states.back(), I);
  if (s) add_old_state(n, t);
  if (s) init_state(old_states.back(), I);
}

void States::update() {
  auto t0 = time();
  for (size_t i = 0; i < old_states.size(); ++i) {
    auto to = old_states[i];
    auto to_name = (std::string)apf::getName(to);
    auto from_name = to_name.erase(to_name.find("_old"), 4);
    auto from = mesh->findField(from_name.c_str());
    apf::copyData(to, from);
  }
  auto t1 = time();
  print(" > states updated in %f seconds", t1 - t0);
}

States* create_states(Disc* d) {
  return new States(d);
}

void destroy_states(States* s) {
  delete s;
}

template void get_scalar(Disc* d, const char* n, apf::MeshEntity* e, ST& v);
template void get_scalar(Disc* d, const char* n, apf::MeshEntity* e, FADT& v);
template void get_tensor(Disc* d, const char* n, apf::MeshEntity* e, minitensor::Tensor<ST>& v);
template void get_tensor(Disc* d, const char* n, apf::MeshEntity* e, minitensor::Tensor<FADT>& v);
template void set_scalar(Disc* d, const char* n, apf::MeshEntity* e, ST const& v);
template void set_scalar(Disc* d, const char* n, apf::MeshEntity* e, FADT const& v);
template void set_tensor(Disc* d, const char* n, apf::MeshEntity* e, minitensor::Tensor<ST> const& v);
template void set_tensor(Disc* d, const char* n, apf::MeshEntity* e, minitensor::Tensor<FADT> const& v);

}
