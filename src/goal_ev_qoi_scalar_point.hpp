#ifndef GOAL_EV_QOI_SCALAR_POINT_HPP
#define GOAL_EV_QOI_SCALAR_POINT_HPP

/** \file goal_ev_qoi_scalar_point.hpp */

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"

/** \cond */
namespace apf {
class Mesh;
class MeshEntity;
}
/** \endcond */

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
class Indexer;
/** \endcond */

/** \brief Compute a point-wise functional for a scalar DOF field.
  * \details This evaluator will compute the value of a scalar DOF field
  * at a given point as \f$ u(x_0) \f$ and appropriately modify the
  * functional derivative vector dJdu with a delta function corresponding
  * to this functional quantity. The functional is given as:
  *
  * \f[
  * J(u) = \int_{\Omega} u \delta(x - x_0) \, \text{d} \Omega.
  * \f]
  *
  * Here \f$ x_0 \f$ is represented by a geometric vertex, such that the
  * mesh generation software will place a mesh vertex that is classified
  * on the geometric vertex at the point \f$ x_0 \f$.
  *
  * dependent fields  | data layout
  * ----------------  | -----------
  * u                 | (Elem, Node)
  *
  * evaluated fields  | data layout
  * ----------------  | -----------
  * pw_u              | (Elem)
  *
  * field descriptions:
  * - u, the scalar DOF \ref goal::Field.
  * - pw_u, the value of the functional integrated over elements. This is
  * zero everywhere, as this evaluator directly modifies the linear algebra
  * data rather than relying on the scatter functional evaluator to
  * populate the fields. */
template <typename EVALT, typename TRAITS>
class QoIScalarPoint : public PHX::EvaluatorWithBaseImpl<TRAITS>,
                       public PHX::EvaluatorDerived<EVALT, TRAITS> {
 public:
  /** \cond */
  typedef typename TRAITS::SetupData SetupData;
  typedef typename TRAITS::PreEvalData PreEvalData;
  typedef typename TRAITS::PostEvalData PostEvalData;
  typedef typename TRAITS::EvalData EvalData;
  typedef typename EVALT::ScalarT ScalarT;
  /** \endcond */

  /** \brief Construct the evaluator. */
  QoIScalarPoint(RCP<Field> u, RCP<Indexer> i, std::string const& set);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

  /** \brief Apply global post-processing for the functional value. */
  void postEvaluate(PostEvalData data);

 private:
  RCP<Field> field;
  RCP<Indexer> indexer;
  std::string set;
  apf::Mesh* mesh;
  apf::MeshEntity* vtx;
  double J;
  PHX::MDField<const ScalarT, Elem, Node> u;
  PHX::MDField<ScalarT, Elem> pw_u;
};

}  // namespace goal

#endif
