#ifndef goal_sol_info_hpp
#define goal_sol_info_hpp

/// @file goal_sol_info.hpp

#include "goal_data_types.hpp"

namespace goal {

/// @cond
class Field;
class Indexer;
class Discretization;
/// @endcond

/// @brief A general linear object container.
struct LinearObj {
  /// @brief The solution increment vector.
  RCP<Vector> du;
  /// @brief The residual vector.
  RCP<Vector> R;
  /// @param The Jacobian matrix.
  RCP<Matrix> dRdu;
  /// @param The functional derivative vector.
  RCP<MultiVector> dJdu;
};

/// @brief A container for parallel solution information.
class SolInfo {

  public:

    /// @brief Construct the solution information object.
    /// @param i The relevant \ref goal::Indexer object.
    /// @param num_qoi The number of QoIs.
    SolInfo(Indexer* i, int num_qoi);

    /// @brief Destroy the solution info object.
    /// @details This will destroy the owned and ghost objects.
    ~SolInfo();

    /// @brief Transfer data from ghost->du to owned->du.
    /// @details This is called with the Tpetra::INSERT directive.
    void gather_du();

    /// @brief Transfer data from ghost->R to owned->R.
    /// @details This is called with the Tpetra::ADD directive.
    void gather_R();

    /// @brief Transfer data from ghost->dRdu to owned->dRdu.
    /// @details This is called with the Tpetra::ADD directive.
    void gather_dRdu();

    /// @brief Transfer data from ghost->dJdu to owned->dJdu.
    /// @details This is called with the Tpetra::ADD directive.
    void gather_dJdu();

    /// @brief Transfer data from owned->du to ghost->du.
    /// @details This is called with the Tpetra::INSERT directive.
    void scatter_du();

    /// @brief Transfer data from owned->R to ghost->R.
    /// @details This is called with the Tpetra::INSERT directive.
    void scatter_R();

    /// @brief Transfer data from owned->dRdu to ghost->dRdu.
    /// @details This is called with the Tpetra::INSERT directive.
    void scatter_dRdu();

    /// @brief Transfer data from owned->dJdu to ghost->dJdu.
    /// @details This is called with the Tpetra::INSERT directive.
    void scatter_dJdu();

    /// @brief The owned linear algebra data.
    LinearObj* owned;

    /// @brief The ghost linear algebra data.
    LinearObj* ghost;

  private:

    RCP<Import> importer;
    RCP<Export> exporter;
};

} // namespace goal

#endif
