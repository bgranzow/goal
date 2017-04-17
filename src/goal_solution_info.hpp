#ifndef GOAL_SOLUTION_INFO_HPP
#define GOAL_SOLUTION_INFO_HPP

/** \file goal_solution_info.hpp */

#include "goal_data_types.hpp"

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
class Indexer;
class Discretization;
/** \endcond */

/** \brief A general linear object container.
  * \details This container generalizes both owned and ghosted linear
  * algebra data structures that are required for both primal and
  * dual solves. */
struct LinearObj {
  /** \brief The DOF solution increment vector. */
  RCP<Vector> du;

  /** \brief The residual vector. */
  RCP<Vector> R;

  /** \brief The Jacobian matrix. */
  RCP<Matrix> dRdu;

  /** \brief The functional derivative vector. */
  RCP<Vector> dJdu;

  /** \brief The dual solution vector. */
  RCP<Vector> z;
};

/** \brief A container for the solution information required for primal
  * and dual solves. */
class SolutionInfo {
 public:
  /** \brief Construct the solution info container.
    * \param indexer The relevant \ref goal::Indexer. */
  SolutionInfo(RCP<Indexer> indexer);

  /** \brief Transfer data from ghost->du to owned->du.
    * \details This is called with the Tpetra::INSERT directive. */
  void gather_du();

  /** \brief Transfer data from ghost->R to owned->R.
    * \details This is called with the Tpetra::ADD directive. */
  void gather_R();

  /** \brief Transfer data from ghost->dRdu to owned->dRdu.
    * \details This is called with the Tpetra::ADD directive. */
  void gather_dRdu();

  /** \brief Transfer data from ghost->dJdu to owned->dJdu.
    * \details This is called with the Tpetra::ADD directive. */
  void gather_dJdu();

  /** \brief Transfer data from ghost->z to owned->z.
    * \details This is called with the Tpetra::ADD directive. */
  void gather_z();

  /** \brief Transfer data from owned->u to ghost->u.
    * \details This is called with the Tpetra::INSERT directive. */
  void scatter_du();

  /** \brief Transfer data from owned->R to ghost->R.
    * \details This is called with the Tpetra::ADD directive. */
  void scatter_R();

  /** \brief Transfer data from owned->dRdu to ghost->dRdu.
    * \details This is called with the Tpetra::ADD directive. */
  void scatter_dRdu();

  /** \brief Transfer data from owned->dRdu to ghost->dRdu.
    * \details This is called with the Tpetra::ADD directive. */
  void scatter_dJdu();

  /** \brief Transfer data from owned->dRdu to ghost->dRdu.
    * \details This is called with the Tpetra::INSERT directive. */
  void scatter_z();

  /** \brief The owned linear algebra data. */
  RCP<LinearObj> owned;

  /** \brief The ghost linear algebra data. */
  RCP<LinearObj> ghost;

 private:
  RCP<Import> importer;
  RCP<Export> exporter;
};

/** \brief Add to fields from the global increment vector.
  * \param f The fields to fill in.
  * \param i The relevant DOF indexer.
  * \param du The relevant solution increment vector. */
void add_to_fields(
    std::vector<RCP<Field> > const& f, RCP<Indexer> i, RCP<Vector> du);

}  /* namespace goal */

#endif
