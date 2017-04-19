#ifndef GOAL_SIZE_FIELD_HPP
#define GOAL_SIZE_FIELD_HPP

/** \file goal_size_field.hpp */

/** \cond */
namespace apf {
class Field;
}
/** \endcond */

namespace goal {

/** \brief Get an isotropic size field targetting an output # of elements.
  * \param e The scalar error field defined at mesh vertices.
  * \param t The target number of elements.
  * \param p The polynomial order of convergence for the FE scheme. */
apf::Field* get_iso_target_size(apf::Field* e, std::size_t t, int p);

} /* namesapce goal */

#endif
