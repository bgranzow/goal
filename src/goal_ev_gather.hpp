#ifndef goal_ev_gather_hpp
#define goal_ev_gather_hpp

/// @file goal_ev_gather.hpp

#include <Phalanx_Evaluator_Derived.hpp>
#include <Phalanx_Evaluator_WithBaseImpl.hpp>

#include "goal_dimension.hpp"
#include "goal_traits.hpp"

namespace goal {

/// @cond
class Field;
class Indexer;
template <typename EVALT, typename TRAITS> class Gather;
/// @endcond

/// @brief Gather from APF fields, Residual specialization.
template <typename TRAITS>
class Gather<goal::Traits::Residual, TRAITS>
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
      public PHX::EvaluatorDerived<goal::Traits::Residual, TRAITS> {

  public:

    /// @cond
    using SetupData = typename TRAITS::SetupData;
    using EvalData = typename TRAITS::EvalData;
    using ScalarT = typename goal::Traits::Residual::ScalarT;
    /// @endcond

    /// @brief Construct the evaluator.
    /// @param i The relevant \ref goal::Indexer.
    /// @param f The \ref goal::Field s to gather data from.
    /// @param t The entity type to operate on.
    Gather(Indexer* i, std::vector<Field*> const& f, int t);

    /// @brief Finalize the field manager registration.
    /// @param d The setup data (void*).
    /// @param fm The phalanx field manager to register.
    void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

    /// @brief Fill in the local multi-dimensional arrays.
    /// @param workset The workset to operate over.
    void evaluateFields(EvalData workset);

  private:

    int num_fields;
    int num_nodes;
    Indexer* indexer;
    std::vector<Field*> fields;
    std::vector<PHX::MDField<ScalarT, Ent, Node> > u;
};

/// @brief Gather from APF fields, Jacobian specialization.
template <typename TRAITS>
class Gather<goal::Traits::Jacobian, TRAITS>
    : public PHX::EvaluatorWithBaseImpl<TRAITS>,
      public PHX::EvaluatorDerived<goal::Traits::Jacobian, TRAITS> {

  public:

    /// @cond
    using SetupData = typename TRAITS::SetupData;
    using EvalData = typename TRAITS::EvalData;
    using ScalarT = typename goal::Traits::Jacobian::ScalarT;
    /// @endcond

    /// @brief Construct the evaluator.
    /// @param i The relevant \ref goal::Indexer.
    /// @param f The \ref goal::Field s to gather data from.
    /// @param t The entity type to operate on.
    Gather(Indexer* i, std::vector<Field*> const& f, int t);

    /// @brief Finalize the field manager registration.
    /// @param d The setup data (void*).
    /// @param fm The phalanx field manager to register.
    void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

    /// @brief Fill in the local multi-dimensional arrays.
    /// @param workset The workset to operate over.
    void evaluateFields(EvalData workset);

  private:

    int num_fields;
    int num_nodes;
    Indexer* indexer;
    std::vector<Field*> fields;
    std::vector<PHX::MDField<ScalarT, Ent, Node> > u;
};

} // end namespace goal

#endif
