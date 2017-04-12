#include "goal_control.hpp"

#include <PCU.h>
#include <RTC_FunctionRTC.hh>
#include <cassert>
#include <cstdarg>
#include <cstdlib>

namespace goal {

static bool is_initd = false;
static bool goal_initd_pcu = false;
static bool goal_initd_mpi = false;
static PG_RuntimeCompiler::Function evaluator(5);

static void init_expression() {
  evaluator.addVar("double", "x");
  evaluator.addVar("double", "y");
  evaluator.addVar("double", "z");
  evaluator.addVar("double", "t");
  evaluator.addVar("double", "val");
}

void initialize(bool init_mpi, bool init_pcu) {
  if (is_initd) return;
  if (init_mpi) {
    MPI_Init(0, 0);
    goal_initd_mpi = true;
  }
  if (init_pcu) {
    PCU_Comm_Init();
    goal_initd_pcu = true;
  }
  init_expression();
  is_initd = true;
}

void finalize() {
  assert(is_initd);
  if (goal_initd_pcu) PCU_Comm_Free();
  if (goal_initd_mpi) MPI_Finalize();
  is_initd = false;
}

void print(const char* message, ...) {
  assert(is_initd);
  if (PCU_Comm_Self()) return void();
  va_list ap;
  va_start(ap, message);
  vfprintf(stdout, message, ap);
  va_end(ap);
  printf("\n");
}

void fail(const char* why, ...) {
  assert(is_initd);
  va_list ap;
  va_start(ap, why);
  vfprintf(stderr, why, ap);
  va_end(ap);
  printf("\n");
  abort();
}

double eval(std::string const& v, double x, double y, double z, double t) {
  std::string value = "val=" + v;
  evaluator.addBody(value);
  evaluator.varValueFill(0, x);
  evaluator.varValueFill(1, y);
  evaluator.varValueFill(2, z);
  evaluator.varValueFill(3, t);
  evaluator.varValueFill(4, 0.0);
  evaluator.execute();
  return evaluator.getValueOfVar("val");
}

double time() { return PCU_Time(); }

}  // namespace goal
