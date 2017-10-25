#include "gmodel.hpp"

int main()
{
  using namespace gmod;
  default_size = 0.5;
  auto s = new_cube(
      Vector{0, 0, 0},
      Vector{1, 0, 0},
      Vector{0, 1, 0},
      Vector{0, 0, 1});
  write_closure_to_geo(s, "cube.geo");
  write_closure_to_dmg(s, "cube.dmg");
}
