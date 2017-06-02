#ifndef goal_assembly_hpp
#define goal_assembly_hpp

/// @file goal_assembly.hpp

namespace goal {

/// @cond
struct Workset;
class Physics;
class SolInfo;
class Discretization;
/// @endcond

/// @brief Assemble volumetric fields.
/// @param ws The \ref goal::Workset to operate on.
/// @param p The \ref goal::Physics to use for assembly.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to opearte on.
template <typename EvalT>
void assemble_volumetric(
    Workset& ws, Physics* p, SolInfo* i, Discretization* d);

/// @brief Assemble neumann fields.
/// @param ws The \ref goal::Workset to operate on.
/// @param p The \ref goal::Physics to use for assembly.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to opearte on.
template <typename EvalT>
void assemble_neumann(
    Workset& ws, Physics* p, SolInfo* i, Discretization* d);

/// @brief Assemble dirichlet fields.
/// @param ws The \ref goal::Workset to operate on.
/// @param p The \ref goal::Physics to use for assembly.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to opearte on.
template <typename EvalT>
void assemble_dirichlet(
    Workset& ws, Physics* p, SolInfo* i, Discretization* d);

/// @brief Compute the forward values of a physical model.
/// @param p The \ref goal::Physics defining the primal model.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to operate on.
/// @param t_now The current evaluation time.
/// @param t_old The previous evaluation time.
void compute_residual(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old);

/// @brief Compute the derivatives of a physical model.
/// @param p The \ref goal::Physics defining the primal model.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to operate on.
/// @param t_now The current evaluation time.
/// @param t_old The previous evaluation time.
void compute_jacobian(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old);

} // end namespace goal

#endif
