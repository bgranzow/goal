#ifndef GOAL_LOG_HPP
#define GOAL_LOG_HPP

#include <vector>

/** \file goal_log.hpp */

namespace goal {

/** \brief A solution statistics logger.
  * \details This class attempts to print convenient statistics about
  * the solution history to make it easy to plot graphs (e.g. convergence
  * graphs). The input arguments dual (was a dual problem solved to
  * perform goal-oriented error estimation?) and exact (is the exact
  * functional value provided to this class?) control the type of
  * output as shown below:
  *
  * `(! dual) && (! exact)`
  * - `time | iter | pDOFs | Ju_h`
  *
  * `(! dual) && (exact)`
  * - `time | iter | pDOFs | Ju_h | J | E`
  *
  * `(dual) && (! exact)`
  * - `time | iter | pDOFs | Ju_h | dDOFs | E_h | B_h`
  *
  * `(dual) && (exact)`
  * - `time | iter | pDOFs | Ju_h | J | E | dDOFs | E_h | B_h | I | IB`
  *
  * Brief descriptions of these values are provided below:
  * - time, The physical evaluation time of the simulation.
  * - iter, The adaptive cycle iteration number.
  * - pDOFs, The number of DOFs for the primal problem.
  * - Ju_h, The FEM approximation to the functional value.
  * - J, The exact functional value.
  * - E, The exact functional error J - Ju_h
  * - dDOFs, The number of DOFs for the dual problem.
  * - E_h, The goal-oriented error estimate for the functional error.
  * - B_h, The estimated upper bound on the functional error.
  * - I, The effectivity of E_h. I = E_h / E.
  * - IB, The effectivity of B_h. IB = B_h / E. */
class Log {
 public:
  /** \brief Construct the log structure.
    * \param dual Was the dual problem used to estimate functional errors?
    * \param exact Is the exact functional value given to this class?
    * \param J The exact functional value. */
  Log(bool dual = false, bool exact = false, double J = 0.0);

  /** \brief Print the current iteration information.
    * \details This prints the most recent iteration that the log class
    * knows about by using the last values set in the public vectors. */
  void print_current();

  /** \brief Print a full summary of solution information.
    * \details This will print information over the entire solution
    * history that has been provided to the log class. */
  void print_summary();

  /** \brief The evaluation time. */
  std::vector<double> time;

  /** \brief The adaptive iteration number. */
  std::vector<int> iter;

  /** \brief The number of DOFs for the primal problem. */
  std::vector<long> pDOFs;

  /** \brief The number of DOFs for the dual problem. */
  std::vector<long> dDOFs;

  /** \brief The approximated functional value. */
  std::vector<double> Ju_h;

  /** \brief The estimated functional error. */
  std::vector<double> E_h;

  /** \brief The estimated upper bound on the functional error. */
  std::vector<double> B_h;

 private:

  double J;
  bool have_dual;
  bool have_exact_J;

  double E;
  double I;
  double IB;

  void print_banner();
  void pre_validate();
  void compute_error(const int i);
  void compute_effectivities(const int i);
  void print_index(const int i);
};

} /* namespace goal */

#endif
