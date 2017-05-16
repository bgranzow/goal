#ifndef goal_ev_basis_hpp
#define goal_ev_basis_hpp

/// @file goal_ev_basis.hpp

#include <Phalanx_Evaluator_Macros.hpp>

#include "goal_dimension.hpp"

namespace goal {

/// @cond
class Field;
/// @endcond

PHX_EVALUATOR_CLASS(Basis)

  public:

    /// @brief Construct the basis evaluator.
    /// @param f A field that uses the basis of interest.
    /// @param type The entity type to operate on.
    Basis(Field* f, int type);
    
  private:

    Field* field;

    int num_nodes;
    int num_ips;
    int num_dims;

    PHX::MDField<double, Ent, IP> wdv;
    PHX::MDField<double, Ent, Node, IP> bf;
    PHX::MDField<double, Ent, Node, IP, Dim> gbf;

PHX_EVALUATOR_CLASS_END

} // end namespace goal

#endif
