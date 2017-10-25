#include "gmodel.hpp"

int main()
{
  using namespace gmod;
  default_size = 0.2;
  auto s = new_square(
      Vector{0, 0, 0},
      Vector{1, 0, 0},
      Vector{0, 1, 0});
  write_closure_to_geo(s, "square.geo");
  write_closure_to_dmg(s, "square.dmg");
}
