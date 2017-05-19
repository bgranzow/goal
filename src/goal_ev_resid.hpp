#ifndef goal_ev_resid_hpp
#define goal_ev_resid_hpp

/// @file goal_ev_resid.hpp

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_data_types.hpp"
#include "goal_dimension.hpp"
#include "goal_traits.hpp"

namespace goal {

using Teuchos::RCP;

/// @cond
class Field;
class Indexer;
template <typename EVALT, typename TRAITS> class Resid;
/// @endcond

/// @brief Assemble the residual vector
template <typename TRAITS>
class Resid<goal::Traits::Residual, TRAITS>
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
      public PHX::EvaluatorDerived<goal::Traits::Residual, TRAITS> {

  public:

    /// @cond
    using SetupData = typename TRAITS::SetupData;
    using PreEvalData = typename TRAITS::PreEvalData;
    using EvalData = typename TRAITS::EvalData;
    using ScalarT = typename goal::Traits::Residual::ScalarT;
    /// @endcond

    /// @brief Construct the resid evaluator.
    /// @param i The relevant \ref goal::Indexer.
    /// @param f The DOF \ref goal::Field s.
    /// @param t The entity type to operate on.
    Resid(Indexer* i, std::vector<Field*> const& f, int t, bool);

    /// @brief Finalize the field manager registration.
    /// @param d The \ref goal::Traits::SetupData.
    /// @param fm The phalanx field manager to register.
    void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

    /// @brief Set the pre evaluation data
    /// @param i The pre-evaluation data.
    void preEvaluate(PreEvalData i);

    /// @brief Add contributions to the residual vector.
    /// @param workset The local workset to operate over.
    void evaluateFields(EvalData workset);

  private:

    int num_fields;
    int num_nodes;

    Indexer* indexer;
    std::vector<Field*> fields;
    SolInfo* info;

    std::vector<PHX::MDField<const ScalarT, Ent, Node> > resids;
};

/// @brief Assemble the residual vector and Jacobian matrix
template <typename TRAITS>
class Resid<goal::Traits::Jacobian, TRAITS>
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
      public PHX::EvaluatorDerived<goal::Traits::Jacobian, TRAITS> {

  public:

    /// @cond
    using SetupData = typename TRAITS::SetupData;
    using PreEvalData = typename TRAITS::PreEvalData;
    using EvalData = typename TRAITS::EvalData;
    using ScalarT = typename goal::Traits::Jacobian::ScalarT;
    /// @endcond

    /// @brief Construct the resid evaluator.
    /// @param i The relevant \ref goal::Indexer.
    /// @param f The DOF \ref goal::Field s.
    /// @param t The entity type to operate on.
    /// @param is True if the adjoint of the Jacobian should be filled.
    Resid(Indexer* i, std::vector<Field*> const& f, int t, bool is);

    /// @brief Finalize the field manager registration.
    /// @param d The \ref goal::Traits::SetupData.
    /// @param fm The phalanx field manager to register.
    void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

    /// @brief Set the pre evaluation data
    /// @param i The pre-evaluation data.
    void preEvaluate(PreEvalData i);

    /// @brief Add contributions to the Jacobian data.
    /// @param workset The local workset to operate over.
    void evaluateFields(EvalData workset);

  private:

    void scatter_primal(RCP<Vector> R, RCP<Matrix> dRdu, EvalData ws);
    void scatter_adjoint(RCP<Vector> R, RCP<Matrix> dRdu, EvalData ws);

    int num_fields;
    int num_nodes;
    int num_dofs;
    bool is_adjoint;

    Indexer* indexer;
    std::vector<Field*> fields;
    SolInfo* info;

    std::vector<PHX::MDField<const ScalarT, Ent, Node> > resids;
};

} // end namespace goal

#endif
