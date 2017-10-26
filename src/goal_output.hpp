#ifndef goal_output_hpp
#define goal_output_hpp

#include <Teuchos_ParameterList.hpp>

namespace goal {

class Disc;

class Output {
  public:
    Output(ParameterList const& p, Disc* d);
    ~Output();
    void write(const double t, const int iter);
  private:
    void write_vtk(const double t);
    Disc* disc;
    ParameterList params;
    int interval;
    bool turn_off;
    std::string name;
    int pos;
    int index;
};

Output* create_output(ParameterList const& p, Disc* d);
void destroy_output(Output* o);

}

#endif
