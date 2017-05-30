#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_data_types.hpp"
#include "goal_discretization.hpp"
#include "goal_states.hpp"

namespace goal {

States::States(Discretization* d, int q) {
  disc = d;
  q_degree = q;
  apf_mesh = disc->get_apf_mesh();
}

States::~States() {
  for (size_t i = 0; i < states.size(); ++i)
    apf::destroyField(states[i]);
  for (size_t i = 0; i < old_states.size(); ++i)
    apf::destroyField(old_states[i]);
}

static apf::Matrix3x3 get_I(Discretization* d) {
  auto dim = d->get_num_dims();
  static apf::Matrix3x3 I1x1(1, 0, 0, 0, 0, 0, 0, 0, 0);
  static apf::Matrix3x3 I2x2(1, 0, 0, 0, 1, 0, 0, 0, 0);
  static apf::Matrix3x3 I3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
  switch (dim) {
    case 1: return I1x1; break;
    case 2: return I2x2; break;
    case 3: return I3x3; break;
    default: fail("unable to set identity");
  }
}

static void set_I(Discretization* d, apf::Field* f, int t) {
  GOAL_ALWAYS_ASSERT(t == apf::MATRIX);
  auto I = get_I(d);
  auto m = d->get_apf_mesh();
  apf::MeshEntity* elem;
  auto elems = m->begin(d->get_num_dims());
  while ((elem = m->iterate(elems))) {
    auto type = m->getType(elem);
    auto nodes = apf::getShape(f)->countNodesOn(type);
    for (int n = 0; n < nodes; ++n)
      apf::setMatrix(f, elem, n, I);
  }
  m->end(elems);
}

static void init_state(Discretization* d, apf::Field* f, int t, bool I) {
  if (I) set_I(d, f, t);
  else apf::zeroField(f);
}

void States::add(const char* n, int t, bool s, bool I) {
  add_state(n, t);
  init_state(disc, states.back(), t, I);
  if (s) {
    add_old_state(n, t);
    init_state(disc, old_states.back(), t, I);
  }
}

void States::add_state(const char* n, int t) {
  auto dim = disc->get_num_dims();
  auto s = apf::getVoronoiShape(dim, q_degree);
  auto f = apf::createField(apf_mesh, n, t, s);
  states.push_back(f);
}

void States::add_old_state(const char* n, int t) {
  auto dim = disc->get_num_dims();
  auto s = apf::getVoronoiShape(dim, q_degree);
  auto on = (std::string)n + "_old";
  auto f = apf::createField(apf_mesh, on.c_str(), t, s);
  old_states.push_back(f);
}

static apf::Field* swap(apf::Field* f, apf::FieldShape* s) {
  auto m = apf::getMesh(f);
  auto n = (std::string)apf::getName(f);
  auto t = apf::getValueType(f);
  auto nf = apf::createField(m, "tmp", t, s);
  apf::projectField(nf, f);
  apf::destroyField(f);
  apf::renameField(nf, n.c_str());
  return nf;
}

void States::project(int q) {
  auto t0 = time();
  q_degree = q;
  auto dim = disc->get_num_dims();
  auto s = apf::getVoronoiShape(dim, q_degree);
  for (size_t i = 0; i < states.size(); ++i)
    states[i] = swap(states[i], s);
  for (size_t i = 0; i < old_states.size(); ++i)
    old_states[i] = swap(old_states[i], s);
  auto t1 = time();
  print(" > state fields projected in %f seconds", t1 - t0);
}

void States::update() {
  for (std::size_t i = 0; i < old_states.size(); ++i) {
    auto to = old_states[i];
    auto to_name = (std::string)apf::getName(to);
    auto from_name = to_name.erase(to_name.find("_old"), 4);
    auto from = apf_mesh->findField(from_name.c_str());
    apf::copyData(to, from);
  }
}

static double get_val(double v) {
  return v;
}

static double get_val(FadType const& v) {
  return v.val();
}

static void zero(apf::Vector3& v) {
  for (int i = 0; i < 3; ++i)
    v[i] = 0.0;
}

static void zero(apf::Matrix3x3& m) {
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j)
    m[i][j] = 0.0;
}

template <typename T>
void States::set_scalar(char const* name, apf::MeshEntity* e, int n,
    T const& v) {
  auto f = apf_mesh->findField(name);
  apf::setScalar(f, e, n, get_val(v));
}

template <typename T>
void States::set_vector(char const* name, apf::MeshEntity* e, int n,
    minitensor::Vector<T> const& v) {
  apf::Vector3 val;
  zero(val);
  for (std::size_t i = 0; i < v.get_dimension(); ++i)
    val[i] = get_val(v(i));
  auto f = apf_mesh->findField(name);
  apf::setVector(f, e, n, val);
}

template <typename T>
void States::set_tensor(char const* name, apf::MeshEntity* e, int n,
    minitensor::Tensor<T> const& v) {
  apf::Matrix3x3 val;
  zero(val);
  for (std::size_t i = 0; i < v.get_dimension(); ++i)
    for (std::size_t j = 0; j < v.get_dimension(); ++j)
      val[i][j] = get_val(v(i, j));
  auto f = apf_mesh->findField(name);
  apf::setMatrix(f, e, n, val);
}

template <typename T>
void States::get_scalar(char const* name, apf::MeshEntity* e, int n,
    T& v) {
  auto f = apf_mesh->findField(name);
  auto val = apf::getScalar(f, e, n);
  v = (T)val;
}

template <typename T>
void States::get_vector(char const* name, apf::MeshEntity* e, int n,
    minitensor::Vector<T>& v) {
  auto f = apf_mesh->findField(name);
  apf::Vector3 val;
  apf::getVector(f, e, n, val);
  for (std::size_t i = 0; i < v.get_dimension(); ++i)
    v(i) = (T)val[i];
}

template <typename T>
void States::get_tensor(char const* name, apf::MeshEntity* e, int n,
    minitensor::Tensor<T>& v) {
  auto f = apf_mesh->findField(name);
  apf::Matrix3x3 val;
  apf::getMatrix(f, e, n, val);
  for (std::size_t i = 0; i < v.get_dimension(); ++i)
    for (std::size_t j = 0; j < v.get_dimension(); ++j)
      v(i, j) = (T)val[i][j];
}

States* create_states(Discretization* d, int q) {
  return new States(d, q);
}

void destroy_states(States* s) {
  delete s;
}

template void States::set_scalar(
    char const* name, apf::MeshEntity* e, int n, double const& v);
template void States::set_scalar(
    char const* name, apf::MeshEntity* e, int n, FadType const& v);
template void States::set_vector(char const* name, apf::MeshEntity* e,
    int n, minitensor::Vector<double> const& v);
template void States::set_vector(char const* name, apf::MeshEntity* e,
    int n, minitensor::Vector<FadType> const& v);
template void States::set_tensor(char const* name, apf::MeshEntity* e,
    int n, minitensor::Tensor<double> const& v);
template void States::set_tensor(char const* name, apf::MeshEntity* e,
    int n, minitensor::Tensor<FadType> const& v);

template void States::get_scalar(
    char const* name, apf::MeshEntity* e, int n, double& v);
template void States::get_scalar(
    char const* name, apf::MeshEntity* e, int n, FadType& v);
template void States::get_vector(
    char const* name, apf::MeshEntity* e, int n, minitensor::Vector<double>& v);
template void States::get_vector(char const* name, apf::MeshEntity* e,
    int n, minitensor::Vector<FadType>& v);
template void States::get_tensor(
    char const* name, apf::MeshEntity* e, int n, minitensor::Tensor<double>& v);
template void States::get_tensor(char const* name, apf::MeshEntity* e,
    int n, minitensor::Tensor<FadType>& v);

}  /* namespace goal */
