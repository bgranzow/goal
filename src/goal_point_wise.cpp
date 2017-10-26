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
  p.set<int>("idx", 0);
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
  apf::Vector3 disp;
  apf::MeshEntity* vtx = 0;
  auto dMdu = s->ghost->dMdu->get1dViewNonConst();
  auto set = params.get<std::string>("point");
  auto idx = params.get<int>("idx");
  auto nodes = this->disc->get_nodes(set);
  auto mesh = this->disc->get_apf_mesh();
  auto u = mesh->findField("u");
  if (nodes.size() > 0) {
    GOAL_DEBUG_ASSERT(nodes.size() == 1);
    vtx = nodes[0].entity;
  }
  if (vtx && mesh->isOwned(vtx)) {
    LO row = this->disc->get_lid(vtx, 0, idx);
    dMdu[row] = 1.0;
    apf::getVector(u, vtx, 0, disp);
    this->qoi_value = disp[idx];
  }
  PCU_Add_Doubles(&(this->qoi_value), 1);
  this->disc = 0;
}

template class PointWise<ST>;
template class PointWise<FADT>;

}
