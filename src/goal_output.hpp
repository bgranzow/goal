#ifndef goal_output_hpp
#define goal_output_hpp

/// @file goal_output.hpp

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Array.hpp>

namespace goal {

/// @cond
class Discretization;
/// @endcond

/// @brief A class to manage output visualization over multiple steps.
/// @details The output object exists to provide an easy interface to
/// manage dumping out visualization output data in VTK file format over
/// multiple time or load steps, including options to write output at
/// specific intervals.
class Output {

  public:

    /// @brief Construct the output object.
    /// @param p A parameter list to describe output behavior.
    /// parameter name  | parameter type
    /// --------------  | --------------
    /// out file        | std::string
    /// interval        | int
    /// turn off        | bool
    /// interpolate     | Teuchos::Array<std::string>
    ///
    /// parameter descriptions:
    /// - out file, the base name to dump output.
    /// - interval, dump output at every ith interval.
    /// - turn off, turn output completely off.
    /// - interpolate, interpolate fields to geometry field.
    /// @param d The relevant \ref goal::Discretization.
    Output(ParameterList const& p, Discretization* d);

    /// @brief Write output at time t.
    void write(const double t);

  private:

    void write_vtk(const double t);

    ParameterList params;
    Discretization* disc;

    int interval;
    bool turn_off;
    Teuchos::Array<std::string> fnames;

    std::string name;
    int pos;
    int index;
};

/// @brief Create an output object.
/// @param p The output object input parameters.
/// @param d The relevant \ref goal::Discretization.
Output* create_output(ParameterList const& p, Discretization* d);

/// @brief Destroy an output object.
/// @param o The \ref goal::Output object to destroy.
void destroy_output(Output* o);

} // end namespace goal

#endif
