#ifndef goal_ev_dual_hpp
#define goal_ev_dual_hpp

/// @file goal_ev_dual.hpp

#include <Phalanx_Evaluator_Macros.hpp>

#include "goal_dimension.hpp"

namespace goal {

/// @cond
class Field;
/// @endcond

PHX_EVALUATOR_CLASS(Dual)

  public:

    /// @brief Construct the dual weight evaluator.
    /// @param z The dual \ref goal::Field on the coarse space V^H.
    /// @param z_fine The dual \ref goal::Field on the fine space V^H.
    /// @param t The entity type to operate on.
    Dual(
        std::vector<Field*> const&  z,
        std::vector<Field*> const& z_fine,
        int t);

  private:

    int num_fields;
    int num_vtx;
    int num_ips;
    int num_dims;

    std::vector<Field*> z;
    std::vector<Field*> z_fine;

    // output
    std::vector<PHX::MDField<ScalarT, Ent, Node, IP> > w;
    std::vector<PHX::MDField<ScalarT, Ent, Node, IP, Dim> > gw;

PHX_EVALUATOR_CLASS_END

} // end namespace goal

#endif
