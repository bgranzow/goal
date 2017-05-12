#ifndef goal_traits_hpp
#define goal_traits_hpp

/// @file goal_traits.hpp

#include <Phalanx_Traits.hpp>
#include <Phalanx_config.hpp>
#include <Sacado.hpp>
#include <Sacado_mpl_find.hpp>
#include <Sacado_mpl_vector.hpp>

namespace goal {

/// @cond
struct Workset;
class SolutionInfo;
/// @endcond

/// @brief The goal evaluation traits structure.
/// @details C++ traits are used to provide 'implementation details' to
/// generic templated code. Phalanx uses this traits structure to determine
/// different evaluation types ato perform specific evaluations of generic
/// templated code. Currently goal implements two evaluation types:
///
/// evaluation type | scalar type
/// --------------- | -----------
/// Residual        | double
/// Jacobian        | \ref goal::FadType
///
/// Additionally, this class defines data that is passed into Phalanx
/// Field Manager evaluation methods:
///
/// data type     | defined as
/// ---------     | ----------
/// SetupData     | void*
/// PreEvalData   | \ref goal::SolutionInfo&
/// EvalData      | \ref goal::Workset&
/// PostEvalData  | \ref goal::SolutionInfo&
///
/// For more information, see the Phalanx User's Guide \cite PhxUserGuide.
struct Traits : public PHX::TraitsBase {

  /// @brief Real type.
  using RealType = double;

  /// @brief Forward automatic differentiaton type.
  using FadType = Sacado::Fad::SLFad<RealType, Goal_FAD_SIZE> FadType;

  /// @brief Residual evaluation type.
  struct Residual{
    /// @brief Residual scalar type.
    using ScalarT = RealType;
  };

  /// @brief Jacobian evaluation type
  struct Jacobian{
    /// @brief Jacobian scalar type.
    using ScalarT = FadType;
  };

  /// @brief A collection of all evaluation types.
  using EvalTypes = Sacado::mpl::vector<Residual, Jacobian>;

  /// @brief Data passed to PHX::FieldManager::postRegistrationSetup.
  using SetupData = void*;

  /// @brief Data passed to PHX::FieldManager::preEvalaute.
  using PreEvalData = SolutionInfo&;

  /// @brief Data passed to PHX::FieldManager::postEvaluate.
  using PostEvalData = SolutionInfo&;
};

} // end namespace goal

/// @cond
namespace PHX {

template <>
struct eval_scalar_types<goal::Traits::Residual> {
  using type = Sacado::mpl::vector<goal::Traits::RealType>;
};

template <>
struct eval_scalar_types<goal::Traits::Jacobian> {
  using type = Sacado::mpl::vector<
    goal::Traits::FadType, goal::Traits::RealType>;
};

} // end namespace PHX

/// @endcond

#endif
