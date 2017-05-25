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

/// @brief Application mode of Dirichlet boundary conditions.
enum DBCMode {PRE, POST};

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

/// @brief Compute the residual of the primal model.
/// @param p The \ref goal::Physics defining the primal model.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to operate on.
/// @param t_now The current evaluation time.
/// @param t_old The previous evaluation time.
/// @param dbc The \ref goal::DBCMode.
void compute_primal_residual(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old, const int dbc = POST);

/// @brief Compute the Jacobian of the primal model.
/// @param p The \ref goal::Physics defining the primal model.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to operate on.
/// @param t_now The current evaluation time.
/// @param t_old The previous evaluation time.
/// @param dbc The \ref goal::DBCMode.
void compute_primal_jacobian(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old, const int dbc = POST);

/// @brief Compute the Jacobian of the dual model.
/// @param p The \ref goal::Physics defining the dual model.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to operate on.
/// @param t_now The current evaluation time.
/// @param t_old The previous evaluation time.
/// @param dbc The \ref goal::DBCMode.
void compute_dual_jacobian(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old, const int dbc = POST);

/// @brief Compute the error residual.
/// @param p The \ref goal::Physics defining the error model.
/// @param i The \ref goal::Indexer used for linear algebra indexing.
/// @param d The \ref goal::Discretization to operate on.
/// @param t_now The current evaluation time.
/// @param t_old The previous evaluation time.
/// @param dbc The \ref goal::DBCMode.
void compute_error_residual(Physics* p, SolInfo* i, Discretization* d,
    const double t_now, const double t_old, const int dbc = POST);

} // end namespace goal

#endif
