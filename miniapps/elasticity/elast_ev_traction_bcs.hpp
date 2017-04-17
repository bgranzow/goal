#ifndef ELAST_EV_TRACTION_BCS_HPP
#define ELAST_EV_TRACTION_BCS_HPP

/** \file elast_ev_traction_bcs.hpp */

#include <Phalanx_Evaluator_WithBaseImpl.hpp>
#include <Phalanx_Evaluator_Derived.hpp>

#include "goal_dimension.hpp"

/** \cond */
namespace goal {
class Indexer;
class Discretization;
class SolutionInfo;
}
/** \endcond */

namespace elast {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;
/** \endcond */

/** \brief An evaluator for traction boundary conditions.
  * \details Modifies the residual vector to appropriately account for
  * contributions of traction defined over side sets. That is, given a user
  * specified traction \f$ t \f$, the right hand side will be modified by
  * the integrals over mesh sides (edges in 2D and faces in 3D) as:
  * \f[
  * ( t, w )_{\Gamma}
  * \f]
  * where \f$ w \f$ denotes the vector-valued weighting function.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * none              | N/A
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * op                | Dummy
  */
template <typename EVALT, typename TRAITS>
class TractionBCs : public PHX::EvaluatorWithBaseImpl<TRAITS>,
                    public PHX::EvaluatorDerived<EVALT, TRAITS> {
 public:
  /** \cond */
  typedef typename TRAITS::SetupData SetupData;
  typedef typename TRAITS::PreEvalData PreEvalData;
  typedef typename TRAITS::PostEvalData PostEvalData;
  typedef typename TRAITS::EvalData EvalData;
  typedef typename EVALT::ScalarT ScalarT;
  using Elem = goal::Elem;
  using Node = goal::Node;
  using Dim = goal::Dim;
  using IP = goal::IP;
  using Dummy = goal::Dummy;
  /** \endcond */

  /** \brief Construct the traction BCs evaluator.
    * \param i The relevant \ref goal::Indexer structure.
    * \param p The traction boundary condition parameterlist, containing
    * entries of type Teuchos::Array<std::string> of the form:
    *  - [ side set name, t_x, t_y, t_z ]. */
  TractionBCs(RCP<goal::Indexer> i, RCP<const ParameterList> p);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Grab the linear algebra data structures.
    * \param info The PreEvalData structure (\ref goal::SolutionInfo). */
  void preEvaluate(PreEvalData info);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  void validate_params();
  void apply_bc(EvalData workset, Teuchos::Array<std::string> const& a);

  RCP<const ParameterList> params;
  RCP<goal::Indexer> indexer;
  RCP<goal::Discretization> disc;
  RCP<goal::SolutionInfo> info;
  int num_dims;
};

}  /* namespace elast */

#endif
