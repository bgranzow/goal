#ifndef goal_control_hpp
#define goal_control_hpp

/// @file goal_control.hpp

#include <string>

namespace goal {

/// @brief Initialize parallel and expression parsing services.
/// @param init_mpi Should a call be made to initialize MPI?
/// @param init_pcu Should a call be made to initialize PCU?
/// @details This method calls initialize routines for MPI and PCU,
/// as well as initializing a real-time string expression parser.
void initialize(bool init_mpi = true, bool init_pcu = true);

/// @brief Finalize the parallel services.
/// @details This method calls finalize routines for MPI and PCU
/// if they were initialized with \ref goal::initialize.
void finalize();

/// @brief Print a printf-style formatted message on rank 0.
void print(const char* msg, ...);

/// @brief Fail the application with an explanation message.
void fail(const char* why, ...) __attribute__((noreturn));

/// @brief Fail the application by assertion, with an error message.
void assert_fail(const char* why, ...) __attribute__((noreturn));

/// @brief Evaluate a string expression.
/// @param v The input mathematical string expression.
/// @param x The x coordinate location.
/// @param y The y coordinate location.
/// @param z The z coordinate location.
/// @param t The current time.
double eval(
    std::string const& v,
    const double x,
    const double y,
    const double z,
    const double t);

/// @brief Get the wall time in seconds.
double time();

} // end namespace goal

/// @brief Always assert a conditional
#define GOAL_ALWAYS_ASSERT(cond)                  \
  do {                                            \
    if (! (cond)) {                               \
      char omsg[2048];                            \
      sprintf(omsg, "%s failed at %s + %d \n",    \
        #cond, __FILE__, __LINE__);               \
      goal::assert_fail(omsg);                    \
    }                                             \
  } while (0)

/// @brief Always assert a conditional with a message
#define GOAL_ALWAYS_ASSERT_VERBOSE(cond, msg)     \
  do {                                            \
    if (! (cond)) {                               \
      char omsg[2048];                            \
      sprintf(omsg, "%s failed at %s + %d \n %s", \
        #cond, __FILE__, __LINE__, msg);          \
      goal::assert_fail(omsg);                    \
    }                                             \
  } while(0)

#ifdef NDEBUG
/// @brief Do nothing - optimized out
#define GOAL_DEBUG_ASSERT(cond)
/// @brief Do nothing - optimized out
#define GOAL_DEBUG_ASSERT_VERBOSE(cond, msg)
#else
/// @brief Assert a conditional
#define GOAL_DEBUG_ASSERT(cond) \
  GOAL_ALWAYS_ASSERT(cond)
/// @brief Assert a conditional
#define GOAL_DEBUG_ASSERT_VERBOSE(cond, msg) \
  GOAL_ALWAYS_ASSERT_VERBOSE(cond, msg)
#endif

#endif
