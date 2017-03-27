#ifndef GOAL_TRAITS_HPP
#define GOAL_TRAITS_HPP

/** \file goal_traits.hpp */

#include <Phalanx_Traits.hpp>
#include <Phalanx_config.hpp>
#include <Sacado.hpp>
#include <Sacado_mpl_find.hpp>
#include <Sacado_mpl_vector.hpp>

namespace goal {

/** \cond */
struct Workset;
class SolutionInfo;
/** \endcond */

/** \brief The Goal Traits structure.
  * \details C++ traits are used to provide 'implementation details' to
  * generic templated code. Phalanx uses this traits structure to determine
  * different evaluation types and the scalar types associated with those
  * evaluation types to peform generic evaluations of templated code.
  * Currently, goal implements two evaluation types:
  *
  * evaluation type | scalar type
  * --------------- | -----------
  * Residual        | double
  * Jacobian        | forward automatic differentiation type
  *
  * Additionally, this structure defines data that is passed into Phalanx
  * evalaute methods:
  *
  * data type     | defined as                | passed to PHX::FieldManager::
  * ---------     | ---------                 | -------------
  * SetupData     | void*                     | postRegistrationSetup
  * PreEvalData   | \ref goal::SolutionInfo&  | preEvaluate
  * EvalData      | \ref goal::Workset&       | evaluateFields
  * PostEvalData  | \ref goal::SolutionInfo&  | postEvaluate
  *
  * For more information, see the Phalanx
  * <a href=https://trilinos.org/docs/dev/packages/phalanx/doc/html/index.html>
  * User's guide</a> */
struct Traits : public PHX::TraitsBase {
  /** \brief Real data type = double. */
  typedef double RealType;

  /** \brief Forward automatic differentiation type. */
  typedef Sacado::Fad::SLFad<RealType, GOAL_FAD_SIZE> FadType;

  /** \brief Residual evaluation type. */
  struct Residual {
    /** \brief Residaul scalar type = double. */
    typedef RealType ScalarT;
  };

  /** \brief Jacobian evaluation type. */
  struct Jacobian {
    /** \brief Jacobian scalar type = FAD type. */
    typedef FadType ScalarT;
  };

  /** \brief A collection of all evaluation types. */
  typedef Sacado::mpl::vector<Residual, Jacobian> EvalTypes;

  /** \brief Data passed into PHX::FieldManager::postRegistrationSetup. */
  typedef void* SetupData;

  /** \brief Data passed into PHX::FieldManager::preEvalaute. */
  typedef SolutionInfo& PreEvalData;

  /** \brief Data passed into PHX::FieldManager::evalData. */
  typedef Workset& EvalData;

  /** \brief Data passed into PHX::FieldManager::postEvaluate. */
  typedef SolutionInfo& PostEvalData;
};

}  // namespace goal

/** \cond */
namespace PHX {

template <>
struct eval_scalar_types<goal::Traits::Residual> {
  typedef Sacado::mpl::vector<goal::Traits::RealType> type;
};

template <>
struct eval_scalar_types<goal::Traits::Jacobian> {
  typedef Sacado::mpl::vector<goal::Traits::FadType, goal::Traits::RealType>
      type;
};
/** \endcond */

}  // namespace PHX

#endif
