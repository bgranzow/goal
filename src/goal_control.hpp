#ifndef goal_control_hpp
#define goal_control_hpp

#include <string>

namespace goal {

void initialize(
    bool init_mpi = true,
    bool init_kokkos = true,
    bool init_pcu = true);

void finalize();

void print(const char* msg, ...);

void fail(const char* why, ...) __attribute__((noreturn));

void assert_fail(const char* why, ...) __attribute__((noreturn));

double eval(
    std::string const& v,
    const double x,
    const double y,
    const double z,
    const double t);

double time();

}

#define GOAL_ALWAYS_ASSERT(cond)                      \
  do {                                                \
    if (! (cond)) {                                   \
      char omsg[2048];                                \
      sprintf(omsg, "%s failed at %s + %d \n",        \
        #cond, __FILE__, __LINE__);                   \
      goal::assert_fail(omsg);                        \
    }                                                 \
  } while (0)

#define GOAL_ALWAYS_ASSERT_VERBOSE(cond, msg)         \
  do {                                                \
    if (! (cond)) {                                   \
      char omsg[2048];                                \
      sprintf(omsg, "%s failed at %s + %d \n %s \n",  \
        #cond, __FILE__, __LINE__, msg);              \
      goal::assert_fail(omsg);                        \
    }                                                 \
  } while(0)

#ifdef NDEBUG
#define GOAL_DEBUG_ASSERT(cond)
#define GOAL_DEBUG_ASSERT_VERBOSE(cond, msg)
#else
#define GOAL_DEBUG_ASSERT(cond) \
  GOAL_ALWAYS_ASSERT(cond)
#define GOAL_DEBUG_ASSERT_VERBOSE(cond, msg) \
  GOAL_ALWAYS_ASSERT_VERBOSE(cond, msg)
#endif

#endif
