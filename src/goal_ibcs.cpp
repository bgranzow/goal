#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_sol_info.hpp"
#include "goal_tbcs.hpp"
#include "goal_vector_weight.hpp"

#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

namespace goal {

using Teuchos::rcp_static_cast;

static ParameterList get_valid_params() {
  ParameterList p;
  p.set<std::string>("side set", "");
  p.set<double>("scale distance", 0.0);
  return p;
}

static void validate_params(ParameterList const& p, SolInfo* s) {
  auto disc = s->get_disc();
  auto ss_name = p.get<std::string>("side set");
  disc->get_sides(ss_name);
  p.validateParameters(get_valid_params(), 0);
}

static void compute_traction(
    double scale,
    apf::Vector3 const& ip,
    apf::Vector3 const& center,
    apf::Vector3& T) {
  T = ip - center;
  auto length = T.getLength();
  T = T/length;
  T = T*scale;
}

static void apply_bc(
    double scale,
    std::string const& ss_name,
    RCP<VectorWeight> w,
    SolInfo* s) {

  auto disc = s->get_disc();
  auto sides = disc->get_sides(ss_name);
  auto num_dims = disc->get_num_dims();
  auto mesh = disc->get_apf_mesh();
  auto shape = mesh->getShape();
  auto R = s->ghost->R->get1dViewNonConst();
  GOAL_DEBUG_ASSERT(num_dims == 3);

  apf::Vector3 xi(0,0,0);
  apf::Vector3 x(0,0,0);
  apf::Vector3 T(0,0,0);
  apf::Vector3 c(0,0,0);

  for (size_t s = 0; s < sides.size(); ++s) {
    auto side = sides[s];
    auto me = apf::createMeshElement(mesh, side);
    auto es = shape->getEntityShape(mesh->getType(side));
    auto num_nodes = es->countNodes();
    w->gather(me);
    w->in_elem(me);
    apf::getIntPoint(me, 1, 0, xi);
    double dv = apf::getDV(me, xi);
    double ipw = apf::getIntWeight(me, 1, 0);
    apf::mapLocalToGlobal(me, xi, x);
    w->at_point(xi, ipw, dv);
    compute_traction(scale, x, c, T);
    for (int n = 0; n < num_nodes; ++n) {
      for (int d = 0; d < num_dims; ++d) {
        LO row = disc->get_lid(side, n, d);
        R[row] -= T[d] * w->val(n, d) * ipw * dv;
      }
    }
    w->out_elem();
    apf::destroyMeshElement(me);
  }

}

void set_ibcs(
    ParameterList const& p,
    RCP<Integrator> w,
    SolInfo* s) {
  if (! p.isParameter("side set")) return;
  validate_params(p, s);
  auto scale = p.get<double>("scale distance");
  auto ss_name = p.get<std::string>("side set");
  auto ww = rcp_static_cast<VectorWeight>(w);
  apply_bc(scale, ss_name, ww, s);
}

}
