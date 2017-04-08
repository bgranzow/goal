#ifndef ELAST_PHYSICS_HPP
#define ELAST_PHYSICS_HPP

/** \file elast_physics.hpp */

#include <goal_physics.hpp>

/** \cond */
namespace goal {
class StateFields;
}
/** \endcond */

/** \brief All elasticity mini-app symbols are contained in this namespace. */
namespace elast {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;
/** \endcond */

/** \brief The elasticity mini-app physics interface.
  * \details This class is responsible for building the necessary data to
  * perform evaluations for the assembly of the primal problem: find
  * \f$ u^H \in V^H \f$ such that:
  * \f[
  * R(u^H, w^H) = 0 \quad \forall w^H \in V^H,
  * \f]
  * the dual problem: find \f$ z^h \in V^h \f$ such that:
  * \f[
  * R'[u^h](z^h, w^h) = J'[u^h](w^h) \quad \forall w^h \in V^h,
  * \f]
  * and the vertex-level error estimates given by:
  * \f[
  * \mathcal{E}_i = | R( u^H, (z^h - z^H) \psi_i) |,
  * \f]
  * where \f$ V^H \f$ and \f$ V^h \f$ denote coarse and fine-scale FEM
  * spaces, and \f$ \psi_i \f$ is a partition of unity realized by linear
  * Lagrange basis functions.
  * For more information, see the \ref Elasticity mini-app documentation. */
class Physics : public goal::Physics {
 private:
  /** \cond */
  using FieldManager = goal::FieldManager;
  using Residual = goal::Traits::Residual;
  using Jacobian = goal::Traits::Jacobian;
  using TRAITS = goal::Traits;
  /** \endcond */

 public:
  /** \brief The physics constructor.
    * \param p A parameter list with the following valid parameters:
    *
    * parameter name  | parameter type
    * --------------  | ---------------
    * dirichlet bcs   | Teuchos::ParameterList
    *
    * parameter descriptions:
    * dirichlet bcs, A parameter list with members of type
    * Teuchos::Array<std::string> of the form:
    * [ dof field index, dof field component, node set name, bc value]
    *
    * \param d The relevant \ref goal::Discretization object. */
  Physics(RCP<const ParameterList> p, RCP<goal::Discretization> d);

  /** \brief Destroy the physics object. */
  ~Physics();

  /** \brief Get the Dirichlet boundary condition parameters. */
  RCP<const ParameterList> get_dbc_params();

 private:
  void set_primal();
  void set_dual();
  void set_error();

  void build_primal_volumetric(FieldManager fm);
  void build_primal_neumann(FieldManager);
  void build_primal_dirichlet(FieldManager fm);
  void build_dual_volumetric(FieldManager fm);
  void build_dual_neumann(FieldManager);
  void build_dual_dirichlet(FieldManager fm);
  void build_error_volumetric(FieldManager fm);
  void build_error_neumann(FieldManager);
  void build_error_dirichlet(FieldManager);

  template <typename EvalT>
  void register_volumetric(FieldManager fm);

  template <typename EvalT>
  void register_neumann(FieldManager fm);

  template <typename EvalT>
  void register_dirichlet(FieldManager fm);

  RCP<const ParameterList> params;
  RCP<goal::StateFields> states;

  bool is_primal;
  bool is_dual;
  bool is_error;

  int q_degree;
  int p_order;
};

}  // namespace elast

#endif
