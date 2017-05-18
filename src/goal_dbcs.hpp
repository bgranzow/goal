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
/// @param f Apply to the fine solution?
/// - true: apply to p->u_fine
/// - false: apply to p->u.
void set_dbc_values(Physics* p, const double t, bool f = false);

/// @brief Apply primal Dirichlet boundary conditions.
/// @param p The relevant \ref goal::Physics structure.
/// @param s The linear algebra data to apply the boundary conditions.
/// @param t The current evaluation time.
/// @param c True if columns should be condensed, preserving symmetry.
template <typename EvalT>
void apply_primal_dbcs(Physics* p, SolInfo* i,
    const double t, bool c = false);

/// @brief Apply dual Dirihclet boundary conditions.
/// @param p The relevant \ref goal::Physics structure.
/// @param s The linear algebra data to apply the boundary conditions.
/// @param t The current evaluation time.
/// @param c True if columns should be condensed, preserving symmetry.
template <typename EvalT>
void apply_dual_dbcs(Physics* p, SolInfo* i,
    const double t, bool c = false);

} // end namespace goal

#endif
