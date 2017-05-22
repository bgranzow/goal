#ifndef goal_dbcs_hpp
#define goal_dbcs_hpp

namespace goal {

/// @cond
class Physics;
class SolInfo;
/// @endcond

/// @brief Prescribe the solution u at DBC nodes.
/// @param p The relevant \ref goal::Physics structure.
/// @param t The current evaluation time.
/// @details This applies DBCs to the fields obtained by
/// indexer->get_field(idx).
void set_dbc_values(Physics* p, const double t);

/// @brief Apply primal Dirichlet boundary conditions.
/// @param p The relevant \ref goal::Physics structure.
/// @param i The linear algebra data to apply the boundary conditions.
/// @param t The current evaluation time.
template <typename EvalT>
void apply_primal_dbcs(Physics* p, SolInfo* i, bool c = false);

/// @brief Apply dual Dirihclet boundary conditions.
/// @param p The relevant \ref goal::Physics structure.
/// @param i The linear algebra data to apply the boundary conditions.
/// @param t The current evaluation time.
template <typename EvalT>
void apply_dual_dbcs(Physics* p, SolInfo* i, bool c = false);

} // end namespace goal

#endif
