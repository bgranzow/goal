#ifndef goal_error_hpp
#define goal_error_hpp

namespace apf {
class Field;
}

namespace goal {

apf::Field* compute_error(apf::Field* ue);
double sum_contribs(apf::Field* ue);

}

#endif
