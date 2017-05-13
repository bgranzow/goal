#ifndef goal_workset_hpp
#define goal_workset_hpp

/// @file goal_workset.hpp

#include <string>
#include <vector>

/// @cond
namespace apf {
class MeshEntity;
}
/// @endcond

namespace goal {

/// @brief The workset structure.
/// @details The workset is a collection of mesh entities defined over
/// a subdomain of the entire finite element domain. This structure is
/// passed to the Phalanx::evaluateFields method to perform local FE
/// analysis.
///
/// For more information, see the Phalanx User's Guide \cite PhxUserGuide.
struct Workset {

  /// @brief The number of mesh entities in this workset.
  int size;

  /// @brief The set name that this workset is associated with.
  std::string set;

  /// @brief The current evaluation time.
  double t_now;

  /// @brief The previous evaluation time.
  double t_old;

  /// @brief A collection of mesh entities, all of the same type.
  std::vector<apf::MeshEntity*> entities;
};

} // end namespace goal

#endif
