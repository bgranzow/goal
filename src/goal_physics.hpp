#ifndef GOAL_PHYSICS_HPP
#define GOAL_PHYSICS_HPP

#include <Phalanx_FieldManager.hpp>

#include "goal_traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ArrayRCP;

/** \brief A field manager templated on \ref goal::Traits. */
typedef RCP<PHX::FieldManager<goal::Traits> > FieldManager;

/** \brief An array of \ref goal::FieldManager 's. */
typedef ArrayRCP<FieldManager> FieldManagers;

/** \cond */
class Field;
class Indexer;
class Discretization;
/** \endcond */

/** \brief The abstract physics class.
  * \details The physics class allows users to implement their owns specific
  * physical models assuming there are three distinct modes of evalaution:
  *
  * - primal -> the forward physical model.
  * - dual -> the dual/adjoint physical model.
  * - error -> the error approximation of the physical model.
  *
  * Each mode of evaluation is responsible for constructing its appropriate
  * Phalanx Field Managers, and registering appropriate Phalanx Evaluators in
  * those Field Managers. Goal assumes that no two states of the physics model
  * exist at the same time. A call to a build* method defines the state of the
  * model.
  *
  * Additionally, these models are assumed to be split into three separate
  * types of evaluations:
  *
  * - volumetric -> evaluations over the interior of the domain (elements)
  * - neumann -> evaluations over the Nuemann boundaries (facets)
  * - dirichlet -> evaluations over the Dirichlet boundaries (nodes)
  *
  * There are separate Phalanx Field Managers associated with each of these
  * three types of evaluations. It is up to the physics model to construct
  * and populate these Field Managers, which are accessed via
  * \ref goal::Physics::get_volumetric, \ref goal::Physics::get_neumann, and
  * \ref goal::Physics::get_dirichlet */
class Physics {
 public:

  /** \brief Construct the physics given a \ref goal::Discretization. */
  Physics(RCP<Discretization> d);

  /** \brief Destroy the physics object.
    * \details This does nothing. */
  virtual ~Physics();

  /** \brief Returns the current indexer. */
  RCP<Indexer> get_indexer() { return indexer; }

  /** \brief Returns the volumetric field managers. */
  FieldManagers get_volumetric() { return volumetric_fms; }

  /** \brief Returns the neumann field manager. */
  FieldManager get_neumann() { return neumann_fm; }

  /** \brief Returns the dirichlet field manager. */
  FieldManager get_dirichlet() { return dirichlet_fm; }

  /** \brief Returns the solution fields on the coarse space. */
  std::vector<RCP<Field> > get_u() { return u; }

  /** \brief Returns the dual solution fields on the coarse space. */
  std::vector<RCP<Field> > get_z() { return z; }

  /** \brief Returns the error fields on mesh vertices (PU). */
  std::vector<RCP<Field> > get_e() { return e; }
  
  /** \brief Returns the solution fields on the fine space. */
  std::vector<RCP<Field> > get_u_fine() { return u_fine; }

  /** \brief Returns the dual solution fields on the fine space. */
  std::vector<RCP<Field> > get_z_fine() { return z_fine; }

  /** \brief Build the indexer corresponding to the coarse space V^H.
    * \details This is built using the fields u. */
  void build_coarse_indexer();

  /** \brief Build the indexer corresponding to the fine space V^h.
    * \details This is built using the fields u_fine. */
  void build_fine_indexer();

  /** \brief Build the indexer corresponding to the partition of unity.
    * \details This is built using the fields e. */
  void build_error_indexer();

  /** \brief Destroy the current indexer. */
  void destroy_indexer();

  /** \brief Build the implemented primal model.
    * \details This requires the protected virtual methods:
    *
    * - build_primal_volumetric
    * - build_primal_neumann
    * - build_primal_dirichlet
    *
    * to be implemented. */
  void build_primal_model();

  /** \brief Build the implemented dual model.
    * \details This requires the protected virtual methods:
    *
    * - build_dual_volumetric
    * - build_dual_neumann
    * - build_dual_dirichlet
    *
    * to be implemented. */
  void build_dual_model();

  /** \brief Build the implemented error model.
    * \details This requires the protected virtual methods:
    *
    * - build_error_volumetric
    * - build_error_neumann
    * - build_error_dirichlet
    *
    * to be implemented. */
  void build_error_model();

  /** \brief Destroy the current model. */
  void destroy_model();

  /** \brief Build the enriched data.
    * \details This will construct the fields:
    *
    * - u_fine
    * - z
    * - z_fine
    * - e */
  void build_enriched_data();

  /** \brief Destroy the enriched data fields. */
  void destroy_enriched_data();

 protected:
  virtual void build_primal_volumetric(FieldManager fm) = 0;
  virtual void build_primal_neumann(FieldManager fm) = 0;
  virtual void build_primal_dirichlet(FieldManager fm) = 0;
  virtual void build_dual_volumetric(FieldManager fm) = 0;
  virtual void build_dual_neumann(FieldManager fm) = 0;
  virtual void build_dual_dirichlet(FieldManager fm) = 0;
  virtual void build_error_volumetric(FieldManager fm) = 0;
  virtual void build_error_neumann(FieldManager fm) = 0;
  virtual void build_error_dirichlet(FieldManager fm) = 0;
  RCP<Discretization> disc;
  RCP<Indexer> indexer;
  FieldManagers volumetric_fms;
  FieldManager dirichlet_fm;
  FieldManager neumann_fm;
  std::vector<RCP<Field> > u;
  std::vector<RCP<Field> > z;
  std::vector<RCP<Field> > e;
  std::vector<RCP<Field> > u_fine;
  std::vector<RCP<Field> > z_fine;
};

/** \brief Set the derivative information for a field manager.
  * \param indexer The relevant \ref goal::Indexer.
  * \param fm The relevant \ref goal::FieldManager.
  * \details This sets the number of total elemental degrees of freedom
  * for the field manager, potentially using the indexer to determine
  * which element block the indexer is currently operating over. */
void set_extended_data_type_dims(RCP<Indexer> indexer, FieldManager fm);

}  // namespace goal

#endif
