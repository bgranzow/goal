#ifndef GOAL_EV_UTILS_HPP
#define GOAL_EV_UTILS_HPP

/** \file goal_ev_utils.hpp */

#include "goal_physics.hpp"

namespace goal {

/** \brief Register a degree of freedom field.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations. */
template <typename EvalT>
void register_dof(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

/** \brief Require the scatter operation for the primal model.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations. */
template <typename EvalT>
void require_primal_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

/** \brief Require the scatter operation for the adjoint model.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param fm The field manager that constructs the physics evaluations. */
template <typename EvalT>
void require_adjoint_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm);

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
  * \details See \ref goal::QoIPNorm. */
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

/** \brief Require the scatter of a point-wise QoI for a scalar DOF field.
  * \param u The relevant DOF \ref goal::Field.
  * \param i The relevant \ref goal::Indexer.
  * \param set The node-set that defines the evaluation point.
  * \param fm The field manager that constructs the physics evaluations.
  * \details See \ref goal::QoIScalarPoint. */
void require_qoi_scalar_point(RCP<Field> u, RCP<Indexer> i,
    std::string const& set, FieldManager fm);

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

}  // namespace goal

#endif
