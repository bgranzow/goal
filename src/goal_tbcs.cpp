#include "goal_control.hpp"
#include "goal_disc.hpp"
#include "goal_sol_info.hpp"
#include "goal_tbcs.hpp"
#include "goal_vector_weight.hpp"

#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>

namespace goal {

using Teuchos::Array;
using Teuchos::getValue;
using Teuchos::rcp_static_cast;

static void validate_params(ParameterList const& p, SolInfo* s) {
  auto disc = s->get_disc();
  auto num_dims = disc->get_num_dims();
  for (auto it = p.begin(); it != p.end(); ++it) {
    auto entry = p.entry(it);
    auto a = getValue<Array<std::string>>(entry);
    auto ss_name = a[0];
    disc->get_sides(ss_name);
    GOAL_DEBUG_ASSERT(a.size() == num_dims + 1);
  }
}

static void apply_bc(
    Array<std::string> const& a,
    SolInfo* s,
    RCP<VectorWeight> w,
    double t) {

  apf::Vector3 x(0,0,0);
  apf::Vector3 xi(0,0,0);
  apf::Vector3 T(0,0,0);

  auto disc = s->get_disc();
  auto ss_name = a[0];
  auto sides = disc->get_sides(ss_name);
  auto num_dims = disc->get_num_dims();
  auto R = s->ghost->R->get1dViewNonConst();
  auto mesh = disc->get_apf_mesh();
  auto shape = mesh->getShape();

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
    for (int d = 0; d < num_dims; ++d)
      T[d] = eval(a[d+1], x[0], x[1], x[2], t);
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

void set_tbcs(
    ParameterList const& p,
    RCP<Integrator> w,
    SolInfo* s,
    double t) {
  validate_params(p, s);
  auto ww = rcp_static_cast<VectorWeight>(w);
  for (auto it = p.begin(); it != p.end(); ++it) {
    auto entry = p.entry(it);
    auto a = getValue<Array<std::string>>(entry);
    apply_bc(a, s, ww, t);
  }
}


}
