#ifndef goal_ev_qoi_hpp
#define goal_ev_qoi_hpp

/// @file goal_ev_qoi.hpp

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"
#include "goal_traits.hpp"

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

/// @cond
class Field;
class Indexer;
template <typename EVALT, typename TRAITS> class QoI;
/// @endcond

/// @brief Do nothing for the residual specialization
template <typename TRAITS>
class QoI<goal::Traits::Residual, TRAITS>
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
      public PHX::EvaluatorDerived<goal::Traits::Residual, TRAITS> {

  public:

    /// @cond
    using SetupData = typename TRAITS::SetupData;
    using PreEvalData = typename TRAITS::PreEvalData;
    using EvalData = typename TRAITS::EvalData;
    using ScalarT = typename goal::Traits::Residual::ScalarT;
    /// @endcond

    /// @brief Construct the QoI evaluator.
    /// @param i The relevant \ref goal::Indexer.
    /// @param f A relevant \ref goal::Field.
    /// @param n The name of the PHX MDField to scatter.
    /// @param type The entity type to operate on.
    QoI(Indexer* i, Field* f, std::string const& n, int type);

    /// @brief Finalize the field manager registration.
    /// @param d The \ref goal::Traits::SetupData.
    /// @param fm The phalanx field manager to register.
    void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

    /// @brief Set the pre evaluation data.
    /// @param i The \ref goal::Traits::PreEvalData.
    void preEvaluate(PreEvalData i);

    /// @brief Do nothing.
    /// @param workset The local workset to operate over.
    void evaluateFields(EvalData workset);
};

/// @brief Scatter the QoI derivative vector.
template <typename TRAITS>
class QoI<goal::Traits::Jacobian, TRAITS>
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
      public PHX::EvaluatorDerived<goal::Traits::Jacobian, TRAITS> {

  public:

    /// @cond
    using SetupData = typename TRAITS::SetupData;
    using PreEvalData = typename TRAITS::PreEvalData;
    using EvalData = typename TRAITS::EvalData;
    using ScalarT = typename goal::Traits::Jacobian::ScalarT;
    /// @endcond

    /// @brief Construct the QoI evaluator.
    /// @param i The relevant \ref goal::Indexer.
    /// @param f A relevant \ref goal::Field.
    /// @param n The name of the PHX MDField to scatter.
    /// @param type The entity type to operate on.
    QoI(Indexer* i, Field* f, std::string const& n, int type);

    /// @brief Finalize the field manager registration.
    /// @param d The \ref goal::Traits::SetupData.
    /// @param fm The phalanx field manager to register.
    void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

    /// @brief Set the pre evaluation data.
    /// @param i The \ref goal::Traits::PreEvalData.
    void preEvaluate(PreEvalData i);

    /// @brief Do nothing.
    /// @param workset The local workset to operate over.
    void evaluateFields(EvalData workset);

  private:

    int num_dofs;

    goal::Indexer* indexer;
    goal::SolInfo* info;

    // input
    PHX::MDField<const ScalarT, Ent> qoi;
};

} // end namespace goal

#endif
