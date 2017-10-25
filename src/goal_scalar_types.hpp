#ifndef goal_scalar_types_hpp
#define goal_scalar_types_hpp

#include <Sacado_Fad_SLFad.hpp>

namespace goal {

using ST = double;
using FADT = Sacado::Fad::SLFad<ST, GOAL_FAD_SIZE>;

}

#endif
