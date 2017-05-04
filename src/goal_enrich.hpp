#ifndef GOAL_ENRICH_HPP
#define GOAL_ENRICH_HPP

#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <apfShape.h>
#include <mthQR.h>
#include <set>
#include <pcu_util.h>
#include <iostream>

namespace enrich {

apf::Field* enrichField(apf::Field* f, int desiredOrder);

}

#endif