#ifndef GOAL_FIELD_TYPES_HPP
#define GOAL_FIELD_TYPES_HPP

/** \file goal_field_types.hpp */

namespace goal {

/** \brief Field types. */
enum FieldType { 
  /** \brief A scalar field type. */
  SCALAR = 0, 
  /** \brief A vector field type. */
  VECTOR = 1,
  /** \brief A tensor field type. */
  TENSOR = 2
};

/** \brief FEM basis types. */
enum BasisType {
  /** \brief Lagrangian basis functions. */
  LAGRANGE = 0,
  /** \brief Hierarchical basis functions. */
  HIERARCHIC = 1
};

}  // namespace goal

#endif
