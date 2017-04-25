#include <apf.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_field.hpp>
#include <goal_indexer.hpp>
#include <goal_solution_info.hpp>
#include <goal_traits.hpp>
#include <goal_workset.hpp>
#include <Phalanx_DataLayout_MDALayout.hpp>

#include "elast_ev_traction_bcs.hpp"

namespace elast {

using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::ParameterEntry;
using Teuchos::getValue;
using Teuchos::Array;

template <typename EVALT, typename TRAITS>
void TractionBCs<EVALT, TRAITS>::validate_params() {
  for (auto it = params->begin(); it != params->end(); ++it) {
    auto entry = params->entry(it);
    auto a = getValue<Array<std::string> >(entry);
    assert(a.size() == num_dims + 1);
    auto set = a[0];
    disc->get_sides(set);
  }
}

template <typename EVALT, typename TRAITS>
TractionBCs<EVALT, TRAITS>::TractionBCs(
    RCP<goal::Indexer> i, RCP<const ParameterList> p) {
  params = p;
  indexer = i;
  assert(indexer->get_num_dof_fields() == 1);
  disc = indexer->get_discretization();
  num_dims = indexer->get_field(0)->get_num_dims();
  validate_params();
  auto name = "Traction BCs";
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename EVALT, typename TRAITS>
void TractionBCs<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  (void)d;
  (void)fm;
}

template <typename EVALT, typename TRAITS>
void TractionBCs<EVALT, TRAITS>::preEvaluate(PreEvalData i) {
  info = rcpFromRef(i);
  assert(info->ghost->R != Teuchos::null);
}

template <typename EVALT, typename TRAITS>
void TractionBCs<EVALT, TRAITS>::apply_bc(
    EvalData workset, Teuchos::Array<std::string> const& a) {
  int idx = 0;
  auto t = workset.t_current;
  auto R = info->ghost->R;
  auto set = a[0];
  auto sides = disc->get_sides(set);
  auto u = indexer->get_field(idx);
  auto q_degree = u->get_q_degree();
  auto basis = u->get_apf_basis();
  auto mesh = indexer->get_apf_mesh();
  apf::Vector3 xi(0, 0, 0);
  apf::Vector3 x(0, 0, 0);
  apf::Vector3 traction(0, 0, 0);
  apf::NewArray<double> BF;
  std::vector<goal::LO> numbers;

  for (size_t i = 0; i < sides.size(); ++i) {
    auto f = sides[i];
    auto me = apf::createMeshElement(mesh, f);
    indexer->get_ghost_lids(f, numbers);
    int num_ips = apf::countIntPoints(me, q_degree);
    auto es = basis->getEntityShape(mesh->getType(f));
    int num_nodes = es->countNodes();
    for (int ip = 0; ip < num_ips; ++ip) {
      apf::getIntPoint(me, q_degree, ip, xi);
      apf::mapLocalToGlobal(me, xi, x);
      apf::getBF(basis, me, xi, BF);
      double w = apf::getIntWeight(me, q_degree, ip);
      double dv = apf::getDV(me, xi);
      for (int i = 0; i < num_dims; ++i) {
        traction[i] = goal::eval(a[i + 1], x[0], x[1], x[2], t);
      }
      int dof = 0;
      for (int node = 0; node < num_nodes; ++node) {
        for (int dim = 0; dim < num_dims; ++dim) {
          goal::LO row = numbers[dof++];
          double val = -BF[node] * traction[dim] * w * dv;
          R->sumIntoLocalValue(row, val);
        }
      }
    }
    apf::destroyMeshElement(me);
  }
}

template <typename EVALT, typename TRAITS>
void TractionBCs<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  for (auto i = params->begin(); i != params->end(); ++i) {
    ParameterEntry const& entry = params->entry(i);
    Array<std::string> a = getValue<Array<std::string> >(entry);
    apply_bc(workset, a);
  }
}

template class TractionBCs<goal::Traits::Residual, goal::Traits>;
template class TractionBCs<goal::Traits::Jacobian, goal::Traits>;

}  /* namespace elast */
