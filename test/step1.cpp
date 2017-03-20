#include <apf.h>
#include <apfBox.h>
#include <goal_control.hpp>

int main() {

  /* [Step1 initialize] */
  goal::initialize();
  /* [Step1 initialize] */

  /** [Step1 finalize] */
  goal::print("Goal: Step 1");
  /** [Step1 finalize] */
  
  /** [Step1 create] */
  auto mesh = apf::makeMdsBox(20, 20, 0, 1.0, 1.0, 0.0, true);
  /** [Step1 createmesh] */

  /** [Step1 verify] */
  mesh->verify();
  /** [Step1 verify] */

  /** [Step1 destroy] */
  apf::destroyMesh(mesh);
  /** [Step1 destroy] */

  /** [Step1 finalize] */
  goal::finalize();
  /** [Step1 finalize] */

  return 0;
}
