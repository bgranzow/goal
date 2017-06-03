#ifndef goal_dbcs_hpp
#define goal_dbcs_hpp

/// @file goal_dbcs.hpp

namespace goal {

/// @cond
class Physics;
class SolInfo;
/// @endcond

/// @brief Prescribe the solution u at DBC nodes.
/// @param p The relevant \ref goal::Physics structure.
/// @param t The current evaluation time.
/// @details This applies DBCs to the fields obtained by
/// p->get_indexer()->get_field(idx).
void set_dbc_values(Physics* p, const double t);

/// @brief Apply primal Dirichlet boundary conditions.
/// @param p The relevant \ref goal::Physics structure.
/// @param i The relevant \ref goal::SolInfo structure.
/// @param t The current evaluation time.
template <typename EvalT>
void apply_dbcs(Physics* p, SolInfo* i, const double t);

} // end namespace goal

#endif
