#include "gmodel.hpp"

int main()
{
  using namespace gmod;
  default_size = 0.1;
  auto s = new_square(
      Vector{0, 0, 0},
      Vector{1, 0, 0},
      Vector{0, 1, 0});
  auto p = new_point2(Vector{0.5,0.5,0});
  embed(s,p);
  write_closure_to_geo(s, "square.geo");
  write_closure_to_dmg(s, "square.dmg");
}
