#include <cmath>
#include <goal_control.hpp>

int main()
{
  goal::initialize();
  goal::print("unit test: control");
  double x = 1.0;
  double y = 2.0;
  double z = 3.0;
  double t = 4.0;
  std::string expr = "sin(x)*cos(y) + 2*exp(t) - z";
  double v1 = goal::eval(expr, x, y, z, t);
  double v2 = sin(1.0)*cos(2.0) + 2.0*exp(t) - z;
  GOAL_ALWAYS_ASSERT(fabs(v1 - v2) < 1.0e-15);
  goal::finalize();
}
