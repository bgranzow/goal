#ifndef goal_physics_hpp
#define goal_physics_hpp

/// @file goal_physics.hpp

#include <Phalanx_FieldManager.hpp>

#include "goal_traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ParameterList;

/// @brief A field manager template on \ref goal::Traits.
using FieldManager = RCP<PHX::FieldManager<goal::Traits> >;
/// @brief An array of \ref goal::FieldManager s.
using FieldManagers = ArrayRCP<FieldManager>;

/// @cond
class Field;
class Indexer;
class Discretization;
/// @endcond

/// @brief The abstract physics class.
/// @details The physics class allows users to implement physics models
/// assuming three distinct evaluation modes:
/// - primal -> the forward physics model.
/// - dual -> the dual/adjoint physical model.
/// - error -> the error localization model.
///
/// Each mode is responsible for constructing its own appropriate
/// Phalanx Field Managers and registering Phalanx Evaluators in these
/// Field Managers. No two physical models may exist at the same time.
///
/// These models are split into three types of evaluations:
/// - volumetric -> evaluations over the interior (elements).
/// - neumann -> evaluations over the Neumann boundaries (sides).
/// - dirichlet -> evaluations over the Dirichlet boundaries (nodes).
class Physics {

  public:

    /// @brief Construct the physics object.
    /// @param d The relevant \ref goal::Discretization.
    Physics(Discretization* d);

    /// @brief Destroy the physics object.
    /// @details This does nothing.
    virtual ~Physics();

    /// @brief Get the dirichlet boundar condition parameters.
    /// @details Each member of this parameter list should be of type:
    /// Teuchos::Array<std::string> with entries of the form:
    /// [ dof name, node set name, bc value].
    virtual ParameterList const& get_dbc_params() = 0;

    /// @brief Build the coarse indexer for the space V^H.
    /// @details This is built using the DOF fields u.
    void build_coarse_indexer();

    /// @brief Build the fine indexer for the space V^h.
    /// @details This is built using the DOF fields u_fine.
    void build_fine_indexer();

    /// @brief Build the error indexer for the partition of unity.
    /// @details This is built using the fields e.
    void build_error_indexer();

    /// @brief Returns the current indexer.
    Indexer* get_indexer() { return indexer; }

    /// @brief Destroy the current indexer.
    void destroy_indexer();

    /// @brief Build the user-implemented primal model.
    void build_primal_model();

    /// @brief Build the user-implemented dual model.
    void build_dual_model();

    /// @brief Build the user-implemented error model.
    void build_error_model();

    /// @brief Destroy the current model.
    void destroy_model();

    /// @brief Returns the volumetric field managers.
    FieldManagers get_volumetric() { return vfms; }

    /// @brief Returns the neumann field managers.
    FieldManagers get_neumann() { return nfms; }

    /// @brief Returns the dirichlet field manager.
    FieldManager get_dirichlet() { return dfm; }

    /// @brief Returns the primal solution on the coarse space V^H.
    std::vector<Field*> const& get_u() { return u; }

    /// @brief Returns the dual solution on the coarse space V^H.
    std::vector<Field*> const& get_z() { return z; }

    /// @brief Returns the error field on the partition of unity.
    std::vector<Field*> const& get_e() { return e; }

    /// @brief Returns the primal solution on the fine space V^h.
    std::vector<Field*> const& get_u_fine() { return u_fine; }

    /// @brief Returns the dual solution on the fine space V^h.
    std::vector<Field*> const& get_z_fine() { return z_fine; }

  protected:

    virtual void build_primal_volumetric(FieldManager) {}
    virtual void build_primal_neumann(FieldManager) {}
    virtual void build_primal_dirichlet(FieldManager) {}

    virtual void build_dual_volumetric(FieldManager) {}
    virtual void build_dual_neumann(FieldManager) {}
    virtual void build_dual_dirichlet(FieldManager) {}

    virtual void build_error_volumetric(FieldManager) {}
    virtual void build_error_neumann(FieldManager) {}
    virtual void build_error_dirichlet(FieldManager) {}

    Discretization* disc;
    Indexer* indexer;

    FieldManagers vfms;
    FieldManagers nfms;
    FieldManager dfm;

    std::vector<Field*> u;
    std::vector<Field*> z;
    std::vector<Field*> e;

    std::vector<Field*> u_fine;
    std::vector<Field*> z_fine;

    int elem_set;
    int side_set;
};

/// @brief Set the derivative information for a field manager.
/// @param i The relevant \ref goal::Indexer.
/// @param fm The relevant \ref goal::FieldManager.
/// @param type The entity type to operate on.
void set_extended_data_type_dims(Indexer* i, FieldManager fm, int type);

} // end namespace goal

#endif
