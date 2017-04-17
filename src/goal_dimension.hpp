#ifndef GOAL_DIMENSION_HPP
#define GOAL_DIMENSION_HPP

/** \file goal_dimension.hpp */

#include <Phalanx_DimTag.hpp>

namespace goal {

PHX_DIM_TAG_DECLARATION(Elem)
PHX_DIM_TAG_DECLARATION(Node)
PHX_DIM_TAG_DECLARATION(Dim)
PHX_DIM_TAG_DECLARATION(IP)
PHX_DIM_TAG_DECLARATION(Dummy)

}  /* namespace goal */

#endif
