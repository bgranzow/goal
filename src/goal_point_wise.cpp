#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <PCU.h>

#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_point_wise.hpp"
#include "goal_sol_info.hpp"

namespace goal {

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("type", "");
  p.set<std::string>("point", "");
  return p;
}

template <typename T>
PointWise<T>::PointWise(ParameterList const& p) {
  params = p;
  params.validateParameters(get_valid_params(), 0);
  this->name = "point wise";
}

template <typename T>
void PointWise<T>::pre_process(SolInfo* s) {
  this->qoi_value = 0.0;
  this->disc = s->get_disc();
}

template <typename T>
void PointWise<T>::post_process(SolInfo* s) {
  apf::MeshEntity* vtx = 0;
  auto dMdu = s->ghost->dMdu;
  auto set = params.get<std::string>("point");
  auto nodes = this->disc->get_nodes(set);
  auto mesh = this->disc->get_apf_mesh();
  auto u = mesh->findField("u");
  if (nodes.size() > 0) {
    GOAL_DEBUG_ASSERT(nodes.size() == 1);
    vtx = nodes[0].entity;
  }
  if (vtx && mesh->isOwned(vtx)) {
    GO row = this->disc->get_lid(vtx, 0, 0);
    dMdu->replaceLocalValue(row, 1.0);
    this->qoi_value = apf::getScalar(u, vtx, 0);
  }
  PCU_Add_Doubles(&(this->qoi_value), 1);
  this->disc = 0;
}

template class PointWise<ST>;
template class PointWise<FADT>;

}
