#ifndef POISSON_PHYSICS_HPP
#define POISSON_PHYSICS_HPP

/** \file poisson_physics.hpp */

#include <goal_physics.hpp>

/** \brief All poisson mini-app symbols are contained in this namespace. */
namespace poisson {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;
/** \endcond */

/** \brief The Poisson mini-app physics interface.
  * \details This class is responsible for building the necessary data to
  * perform evaluations for the assembly of the primal problem: find
  * \f$ u^H \in V^H \f$ such that:
  * \f[
  * (\nabla w^H, \nabla u^H) = (w^H, f) \quad \forall w^H \in V^H,
  * \f]
  * the dual problem: find \f$ z^h \in V^h \f$ such that:
  * \f[
  * (\nabla z^h, \nabla w^h) = (w^h, j) \quad \forall w^h \in V^h,
  * \f]
  * and the vertex-level error estimates given by
  * \f[
  * \mathcal{E}_i = | ( \nabla ( (z^h - z^H) \psi_i), u^H) -
  * ((z^h - z^H) \psi_i, f) |
  * \f]
  * where \f$ V^H \f$ and \f$ V^h \f$ denote coarse and fine-scale FEM
  * spaces, and \f$ \psi_i \f$ is a partition of unity realized by linear
  * Lagrange basis functions.
  * For more information, see the \ref poisson mini-app documentation. */
class Physics : public goal::Physics {
 private:
   /** \cond */
   using FieldManager = goal::FieldManager;
   using Residual = goal::Traits::Residual;
   using Jacobian = goal::Traits::Jacobian;
   /** \endcond */

 public:
  /** \brief The physics constructor.
    * \param p A parameter list with the following valid parameters:
    *
    * parameter name    | parameter type
    * ---------------   | --------------
    * forcing function  | std::string
    * functional type   | std::string
    * point set         | std::string
    * dirichlet bcs     | Teuchos::ParameterList
    *
    * parameter descriptions:
    * - forcing function, A mathematical string expression for the right hand
    * side of the primal problem.
    * - functional type, The type of functional quantity of interest. Either
    * "point-wise" or "solution average".
    * - point set, If the functional type "point-wise" is chosen, then this
    * parameter points to an appropriate node set for the point-wise geometric
    * vertex.
    * dirichlet bcs, A parameter list with members of type
    * Teuchos::Array<std::string> of the form:
    * [ dof field index, dof field component, node set name, bc value ]
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
  void build_primal_neumann(FieldManager) {}
  void build_primal_dirichlet(FieldManager fm);
  void build_dual_volumetric(FieldManager fm);
  void build_dual_neumann(FieldManager) {}
  void build_dual_dirichlet(FieldManager fm);
  void build_error_volumetric(FieldManager fm);
  void build_error_neumann(FieldManager) {}
  void build_error_dirichlet(FieldManager) {}

  template <typename EvalT>
  void register_volumetric(FieldManager fm);

  template <typename EvalT>
  void register_dirichlet(FieldManager fm);

  RCP<const ParameterList> params;
  std::string ff;
  std::string functional_type;
  std::string set;

  bool is_primal;
  bool is_dual;
  bool is_error;
};

}  // namespace poisson

#endif
