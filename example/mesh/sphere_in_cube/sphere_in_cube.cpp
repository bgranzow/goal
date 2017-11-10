#include "gmodel.hpp"

using namespace gmod;

Vector vec(double x, double y, double z) {
  return Vector{x,y,z};
}

int main() {
  default_size = 0.2;
  auto outer = new_cube(
      vec(-1,-1,-1),
      vec(2,0,0),
      vec(0,2,0),
      vec(0,0,2));
  auto inner = new_sphere(
      vec(0,0,0),
      vec(1,0,0),
      vec(0,0.5,0));
  auto sphere = new_volume2(inner);
  insert_into(outer, sphere);
  write_closure_to_geo(outer, "sphere_in_cube.geo");
  write_closure_to_dmg(outer, "sphere_in_cube.dmg");
}
