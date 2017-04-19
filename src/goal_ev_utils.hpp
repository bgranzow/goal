#ifndef GOAL_EV_UTILS_HPP
#define GOAL_EV_UTILS_HPP

/** \file goal_ev_utils.hpp */

#include "goal_physics.hpp"

namespace goal {

/** \brief Register a degree of freedom field.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See the evaluators:
  * - \ref goal::GatherScalar<goal::Traits::Residual, TRAITS>
  * - \ref goal::GatherScalar<goal::Traits::Jacobian, TRAITS>
  * - \ref goal::ScalarShape
  * - \ref goal::InterpolateScalar
  * - \ref goal::GatherVector<goal::Traits::Residual, TRAITS>
  * - \ref goal::GatherVector<goal::Traits::Jacobian, TRAITS>
  * - \ref goal::VectorShape
  * - \ref goal::InterpolateVector */
template <typename EvalT>
void register_dof(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

/** \brief Register a dual weighting field.
  * \param z The coarse dual \ref goal::Field on \f$ V^H \f$.
  * \param z_fine The fine dual \ref goal::Field on \f$ V^h \f$.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See the evaluators:
  * - \ref goal::DualScalarWeight
  * - \ref goal::DualVectorWeight */
template <typename EvalT>
void register_dual(RCP<Field> z, RCP<Field> z_fine, FieldManager fm);

/** \brief Require the scatter operation for the primal model.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See the evaluators:
  * - \ref goal::ScatterScalar<goal::Traits::Residual, TRAITS>
  * - \ref goal::ScatterScalar<goal::Traits::Jacobian, TRAITS>
  * - \ref goal::ScatterVector<goal::Traits::Residual, TRAITS>
  * - \ref goal::ScatterVector<goal::Traits::Jacobian, TRAITS> */
template <typename EvalT>
void require_primal_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

/** \brief Require the scatter operation for the adjoint model.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See the evaluators:
  * - \ref goal::ScatterScalar<goal::Traits::Jacobian, TRAITS>
  * - \ref goal::ScatterVector<goal::Traits::Jacobian, TRAITS> */
template <typename EvalT>
void require_adjoint_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

/** \brief Require primal Dirichlet boundary conditions.
  * \param p The relevant DBC parameterlist.
  * \param i The relevant \ref goal::Indexer.
  * \param condense True if DBC columns should be condensed from the Jacobian.
  * \param adj True if the dual boundary conditions should be applied.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See the evaluator:
  * - \ref goal::DirichletBCs<goal::Traits::Residual, TRAITS>
  * - \ref goal::DirichletBCs<goal::Traits::Jacobian, TRAITS> */
template <typename EvalT>
void require_dbc(RCP<const ParameterList> p, RCP<Indexer> i,
    bool adj, bool condense, FieldManager fm);

/** \brief Require the scatter of the KS QoI.
  * \param p A QoI parameter list with the following parameters:
  *
  * parameter name  | parameter type
  * --------------  | --------------
  * name            | std::string
  * scalar name     | std::string
  * p               | int
  * m               | int
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See \ref goal::QoIKS. */
void require_qoi_ks(RCP<const ParameterList> p, RCP<Field> u,
    RCP<Indexer> i, FieldManager fm);

/** \brief Require the scatter of the p-norm QoI.
  * \param p A QoI parameter list with the following parameters:
  *
  * parameter name  | parameter type
  * --------------  | --------------
  * name            | std::string
  * scalar name     | std::string
  * p               | int
  * m               | int
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See \ref goal::QoIPNorm. */
void require_qoi_pnorm(RCP<const ParameterList> p, RCP<Field> u,
    RCP<Indexer> i, FieldManager fm);

/** \brief Require the scatter of a QoI from above.
  * \param p A QoI parameter list. The parameter "name" will
  * determine which QoI is required:
  * - "ks", see \ref goal::QoIKS.
  * - "pnorm", see \ref goal::QoIPNorm.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations. */
void require_qoi(RCP<const ParameterList> p, RCP<Field> u, RCP<Indexer> i,
    FieldManager fm);

/** \brief Require the scatter of a point-wise QoI for a scalar DOF field.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param set The node-set that defines the evaluation point.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See \ref goal::QoIScalarPoint. */
void require_qoi_scalar_point(RCP<Field> u, RCP<Indexer> i,
    std::string const& set, FieldManager fm);

/** \brief Require the scatter of a point-wise QoI for a vector DOF field.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param c The component of the DOF \ref goal::Field.
  * \param set The node-set that defines the evaluation point.
  * \param fm The field manager that constructs the physics evalautions.
  * \details See \ref goal::QoIVectorPoint. */
void require_qoi_vector_point(RCP<Field> u, RCP<Indexer> i, int c,
    std::string const& set, FieldManager fm);

/** \brief Require the fill of the error field via dual weighted residual.
  * \param u The relevant DOF \ref goal::Field.
  * \param e The relevant error \ref goal::Field to fill.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See the evaluators:
  * - \ref goal::ScalarError
  * - \ref goal::VectorError */
void require_error(RCP<Field> u, RCP<Field> e, RCP<Indexer> i,
    FieldManager fm);

}  /* namespace goal */

#endif
