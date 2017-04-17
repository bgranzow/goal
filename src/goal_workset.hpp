#ifndef GOAL_WORKSET_HPP
#define GOAL_WORKSET_HPP

/** \file goal_workset.hpp */

#include <string>
#include <vector>

/** \cond */
namespace apf {
class MeshEntity;
}
/** \endcond */

namespace goal {

/** \brief The workset structure.
  * \details The workset is a collection of information defined over a
  * subdomain of the entire finite element domain that is passed to the
  * Phalanx::evaluateFields method, to perform local FE analysis. For more
  * information see the Phalanx
  * <a href=https://trilinos.org/docs/dev/packages/phalanx/doc/html/index.html>
  * User's guide</a> */
struct Workset {
  /** \brief The number of mesh entities in this workset. */
  int size;

  /** \brief The element block name that this workset is associated with. */
  std::string block;

  /** \brief The current evaluation time. */
  double t_current;

  /** \brief The previous evaluation time. */
  double t_previous;

  /** \brief A collection of mesh entities associated with this workset. */
  std::vector<apf::MeshEntity*> entities;
};

}  /* namespace goal */

#endif
