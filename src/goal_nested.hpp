#ifndef goal_nested_hpp
#define goal_nested_hpp

#include "goal_disc.hpp"

namespace apf {
class MeshTag;
}

namespace goal {

enum Mode { UNIFORM, LONG, SINGLE };

using EntArray = std::vector<apf::MeshEntity*>;

class Nested : public Disc {
  public:
    Nested(Disc* d, int mode);
    ~Nested();
  private:
    void number_elems();
    void copy_mesh();
    int mode;
    apf::Mesh2* base_mesh;
    apf::MeshTag* old_vtx_tag;
    apf::MeshTag* new_vtx_tag;
    EntArray base_elems;
    EntArray old_vertices;
};

Nested* create_nested(Disc* d, const int mode);
void destroy_nested(Nested* n);

}

#endif
