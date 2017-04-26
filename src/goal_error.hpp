#ifndef GOAL_ERROR_HPP
#define GOAL_ERROR_HPP

/** \file goal_error.hpp */

#include <Teuchos_RCP.hpp>

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
/** \endcond */

/** \brief Sum all loaclized error contributions.
  * \param e Localized error fields at mesh vertices.
  * \returns The value:
  * \f[
  * \sum_{f=0}^{nf} \sum_{i=0}^{ne} \sum_{c=0}^{nc} \left[ e^f_i \right]_c
  * \f]
  * where \f$ nf \f$ is the number of error fields and \f$ ni \f$ is the
  * number of entities that the field has been localized to, and
  * \f$ nc \f$ is the number of components of the field. */
double sum_contributions(std::vector<RCP<Field> > const& e);

/** \brief Sum absolute values of localized error contributions.
  * \param e Localized error fields at mesh vertices.
  * \returns The value:
  * \f[
  * \sum_{f=0}^{nf} \sum_{i=0}^{ne} | \sum_{c=0}^{nc} \left[ e^f_i \right]_c |
  * \f]
  * where \f$ nf \f$ is the number of error fields and \f$ ni \f$ is the
  * number of entities that the field has been localized to, and
  * \f$ nc \f$ is the number of components of the field. */
double approx_upper_bound(std::vector<RCP<Field> > const& e);

/** \brief Sum absolute values of error contributions at the entity level.
  * \param e Localized error fields corresponding to DOF fields.
  * \returns A scalar field at mesh vertices used to compute a size field. */
apf::Field* compute_indicators(std::vector<RCP<Field> > const& e);

} /* namespace goal */

#endif
