#ifndef goal_states_hpp
#define goal_states_hpp

#include <MiniTensor.h>

namespace apf {
class Field;
class Mesh;
}

namespace goal {

template <typename T>
void get_scalar(Disc* d, const char* n, apf::MeshEntity* e, T& v);

template <typename T>
void get_tensor(Disc* d, const char* n, apf::MeshEntity* e, minitensor::Tensor<T>& v);

template <typename T>
void set_scalar(Disc* d, const char* n, apf::MeshEntity* e, T const& v);

template <typename T>
void set_tensor(Disc* d, const char* n, apf::MeshEntity* e, minitensor::Tensor<T> const& v);

class States {
  public:
    States(Disc* d);
    ~States();
    Disc* get_disc() { return disc; }
    void add(const char* n, int type, bool save_old = false, bool identitize = false);
    void update();
  private:
    void add_state(const char* n, int type);
    void add_old_state(const char* n, int type);
    int num_dims;
    Disc* disc;
    apf::Mesh* mesh;
    std::vector<apf::Field*> states;
    std::vector<apf::Field*> old_states;
};

States* create_states(Disc* d);
void destroy_states(States* s);

}

#endif
