#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

#include "goal_control.hpp"
#include "goal_data_types.hpp"
#include "goal_discretization.hpp"
#include "goal_state_fields.hpp"

namespace goal {

const apf::ValueType apf_types[3] = {apf::SCALAR, apf::VECTOR, apf::MATRIX};

StateFields::StateFields(RCP<Discretization> d, int q) {
  disc = d;
  q_degree = q;
  apf_mesh = disc->get_apf_mesh();
}

static apf::Matrix3x3 get_identity(RCP<Discretization> d) {
  static apf::Matrix3x3 I2x2(1, 0, 0, 0, 1, 0, 0, 0, 0);
  static apf::Matrix3x3 I3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
  if (d->get_num_dims() == 2)
    return I2x2;
  else
    return I3x3;
}

static void set_identity(RCP<Discretization> d, apf::Field* f, int t) {
  assert(t == TENSOR);
  auto I = get_identity(d);
  auto am = d->get_apf_mesh();
  apf::MeshEntity* elem;
  auto it = am->begin(d->get_num_dims());
  while ((elem = am->iterate(it))) {
    auto type = am->getType(elem);
    auto nodes = apf::getShape(f)->countNodesOn(type);
    for (int n = 0; n < nodes; ++n) apf::setMatrix(f, elem, n, I);
  }
  am->end(it);
}

void StateFields::add(const char* n, int t, bool save, bool I) {
  auto dim = disc->get_num_dims();
  auto s = apf::getVoronoiShape(dim, q_degree);
  auto f = apf::createField(apf_mesh, n, apf_types[t], s);
  if (!I)
    apf::zeroField(f);
  else
    set_identity(disc, f, t);
  states.push_back(f);
  if (save) {
    auto oname = std::string(n) + "_old";
    auto g = apf::createField(apf_mesh, oname.c_str(), apf_types[t], s);
    if (!I)
      apf::zeroField(g);
    else
      set_identity(disc, g, t);
    states.push_back(g);
    old_states.push_back(g);
  }
}

static double get_val(double v) { return v; }

static double get_val(FadType const& v) { return v.val(); }

static void zero(apf::Vector3& v) {
  for (int i = 0; i < 3; ++i) v[i] = 0.0;
}

static void zero(apf::Matrix3x3& m) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) m[i][j] = 0.0;
}

template <typename T>
void StateFields::set_scalar(
    char const* name, apf::MeshEntity* e, int n, T const& v) {
  auto f = apf_mesh->findField(name);
  apf::setScalar(f, e, n, get_val(v));
}

template <typename T>
void StateFields::set_vector(char const* name, apf::MeshEntity* e, int n,
    minitensor::Vector<T> const& v) {
  apf::Vector3 val;
  zero(val);
  for (std::size_t i = 0; i < v.get_dimension(); ++i) val[i] = get_val(v(i));
  auto f = apf_mesh->findField(name);
  apf::setVector(f, e, n, val);
}

template <typename T>
void StateFields::set_tensor(char const* name, apf::MeshEntity* e, int n,
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
void StateFields::get_scalar(
    char const* name, apf::MeshEntity* e, int n, T& v) {
  auto f = apf_mesh->findField(name);
  auto val = apf::getScalar(f, e, n);
  v = (T)val;
}

template <typename T>
void StateFields::get_vector(
    char const* name, apf::MeshEntity* e, int n, minitensor::Vector<T>& v) {
  auto f = apf_mesh->findField(name);
  apf::Vector3 val;
  apf::getVector(f, e, n, val);
  for (std::size_t i = 0; i < v.get_dimension(); ++i) v(i) = (T)val[i];
}

template <typename T>
void StateFields::get_tensor(
    char const* name, apf::MeshEntity* e, int n, minitensor::Tensor<T>& v) {
  auto f = apf_mesh->findField(name);
  apf::Matrix3x3 val;
  apf::getMatrix(f, e, n, val);
  for (std::size_t i = 0; i < v.get_dimension(); ++i)
    for (std::size_t j = 0; j < v.get_dimension(); ++j) v(i, j) = (T)val[i][j];
}

static apf::Field* swap_field(apf::Field* f, apf::FieldShape* s) {
  auto m = apf::getMesh(f);
  auto name = (std::string)(apf::getName(f));
  auto vt = apf::getValueType(f);
  auto nf = apf::createField(m, "tmp", vt, s);
  apf::zeroField(nf);
  apf::projectField(nf, f);
  apf::destroyField(f);
  apf::renameField(nf, name.c_str());
  return nf;
}

void StateFields::project(int q) {
  auto t0 = time();
  auto dim = disc->get_num_dims();
  q_degree = q;
  auto s = apf::getVoronoiShape(dim, q_degree);
  for (std::size_t i = 0; i < states.size(); ++i)
    states[i] = swap_field(states[i], s);
  auto t1 = time();
  print(" > state fields projected in %f seconds", t1 - t0);
}

void StateFields::update() {
  for (std::size_t i = 0; i < old_states.size(); ++i) {
    auto to = old_states[i];
    auto to_name = (std::string)apf::getName(to);
    auto from_name = to_name.erase(to_name.find("_old"), 4);
    auto from = apf_mesh->findField(from_name.c_str());
    apf::copyData(to, from);
  }
}

template void StateFields::set_scalar(
    char const* name, apf::MeshEntity* e, int n, double const& v);
template void StateFields::set_scalar(
    char const* name, apf::MeshEntity* e, int n, FadType const& v);
template void StateFields::set_vector(char const* name, apf::MeshEntity* e,
    int n, minitensor::Vector<double> const& v);
template void StateFields::set_vector(char const* name, apf::MeshEntity* e,
    int n, minitensor::Vector<FadType> const& v);
template void StateFields::set_tensor(char const* name, apf::MeshEntity* e,
    int n, minitensor::Tensor<double> const& v);
template void StateFields::set_tensor(char const* name, apf::MeshEntity* e,
    int n, minitensor::Tensor<FadType> const& v);

template void StateFields::get_scalar(
    char const* name, apf::MeshEntity* e, int n, double& v);
template void StateFields::get_scalar(
    char const* name, apf::MeshEntity* e, int n, FadType& v);
template void StateFields::get_vector(
    char const* name, apf::MeshEntity* e, int n, minitensor::Vector<double>& v);
template void StateFields::get_vector(char const* name, apf::MeshEntity* e,
    int n, minitensor::Vector<FadType>& v);
template void StateFields::get_tensor(
    char const* name, apf::MeshEntity* e, int n, minitensor::Tensor<double>& v);
template void StateFields::get_tensor(char const* name, apf::MeshEntity* e,
    int n, minitensor::Tensor<FadType>& v);

}  // namespace goal
