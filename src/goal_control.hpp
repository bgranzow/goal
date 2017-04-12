#ifndef GOAL_CONTROL_HPP
#define GOAL_CONTROL_HPP

/** \file goal_control.hpp */

#include <string>

namespace goal {

/** \brief Initialize parallel services.
  * \param init_mpi Should a call be made to initialize MPI?
  * \param init_pcu Should a call be made to initialize PCU?
  * \details This method calls initialize routines for MPI and PCU, as well
  * as initializing the real-time expression string parser. */
void initialize(bool init_mpi = true, bool init_pcu = true);

/** \brief Finalize parallel services.
  * \details This method calls finalize routines for MPI and PCU. */
void finalize();

/** \brief Print a printf-style formatted message on process 0. */
void print(const char* message, ...);

/** \brief Fail the application with an explanation message. */
void fail(const char* why, ...) __attribute__((noreturn));

/** \brief Evaluate a string expression.
  * \param v The input mathematical string expression.
  * \param x The x coordinate location.
  * \param y The y coordinate location.
  * \param z The z coordinate location.
  * \param t The current time. */
double eval(std::string const& v, double x, double y, double z, double t);

/** \brief Get the wall time in seconds. */
double time();

}  // namespace goal

#endif
