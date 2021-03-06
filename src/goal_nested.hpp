#ifndef goal_nested_hpp
#define goal_nested_hpp

#include "goal_disc.hpp"

namespace apf {
class Field;
class MeshTag;
}

namespace goal {

enum Mode { UNIFORM, LONG, SINGLE };

using EntArray = std::vector<apf::MeshEntity*>;

class Nested : public Disc {
  public:
    Nested(Disc* d, int mode);
    ~Nested();
    void set_fine(RCP<VectorT> z, apf::Field* u, apf::Field* p);
    void set_coarse(apf::Field* u, apf::Field* p);
    apf::Field* set_error(apf::Field* nested_err);
  private:
    void number_elems();
    void copy_mesh();
    void tag_old_verts();
    void refine_uniform();
    void refine_mesh();
    void refine_long();
    void refine_single();
    void store_old_verts();
    void subtract_adj(apf::Field* zu, apf::Field* zp);
    void zero_adj(apf::Field* zu, apf::Field* zp);
    int mode;
    double ratio;
    apf::Mesh2* base_mesh;
    apf::MeshTag* old_vtx_tag;
    apf::MeshTag* new_vtx_tag;
    EntArray base_elems;
    EntArray old_vertices;
};

Nested* create_nested(Disc* d, int mode);
void destroy_nested(Nested* n);

}

#endif
