#ifndef GOAL_FIELD_TYPES_HPP
#define GOAL_FIELD_TYPES_HPP

/** \file goal_field_types.hpp */

namespace goal {

/** \brief Field types. */
enum FieldType { SCALAR = 0, VECTOR = 1, TENSOR = 2 };

/** \brief FEM basis types. */
enum BasisType { LAGRANGE = 0, HIERARCHIC = 1 };

}  // namespace goal

#endif
