#ifndef goal_dimension_hpp
#define goal_dimension_hpp

/// @file goal_dimension.hpp

#include <Phalanx_DimTag.hpp>

namespace goal {

PHX_DIM_TAG_DECLARATION(Ent)
PHX_DIM_TAG_DECLARATION(Node)
PHX_DIM_TAG_DECLARATION(Dim)
PHX_DIM_TAG_DECLARATION(IP)
PHX_DIM_TAG_DECLARATION(Dummy)

} // end namespace goal

#endif
