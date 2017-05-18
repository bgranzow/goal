#ifndef goal_ev_interpolate_hpp
#define goal_ev_interpolate_hpp

/// @file goal_ev_interpolate.hpp

#include <Phalanx_Evaluator_Macros.hpp>

#include "goal_dimension.hpp"

namespace goal {

/// @cond
class Field;
/// @endcond

PHX_EVALUATOR_CLASS(Interpolate)

  public:

    /// @brief Construct the interpolation evaluator.
    /// @param f The \ref goal::Field s to interpolate.
    /// @param t The entity type to interpolate onto ips.
    Interpolate(std::vector<Field*> const& f, int t);

  private:

    int num_fields;
    int num_nodes;
    int num_ips;
    int num_dims;
    
    // input
    PHX::MDField<const double, Ent, Node, IP> bf;
    PHX::MDField<const double, Ent, Node, IP, Dim> gbf;
    std::vector<PHX::MDField<const ScalarT, Ent, Node> > nodal;

    // output
    std::vector<PHX::MDField<ScalarT, Ent, IP> > u;
    std::vector<PHX::MDField<ScalarT, Ent, IP, Dim> > gu;

PHX_EVALUATOR_CLASS_END

} // end namespace goal

#endif
