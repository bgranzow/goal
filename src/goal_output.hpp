#ifndef GOAL_OUTPUT_HPP
#define GOAL_OUTPUT_HPP

/** \file goal_output.hpp */

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

/** \cond */
namespace Teuchos {
class ParameterList;
};
/** \endcond */

namespace goal {

/** \cond */
using Teuchos::RCP;

class Discretization;
/** \endcond */

/** \brief A class to manage visualization output over multiple steps.
  * \details The output object exists to provide an easy interface to manage
  * dumping out visualization output data in VTK file format over multiple
  * time or load steps, including options to write output only at specific
  * intervals and to turn off output completely. */
class Output {
 public:
  /** brief The output constructor.
    * \param p A parameterlist to describe the output behavior
    *
    * Parameter Name  | Parameter Type
    * --------------- | --------------
    * out file        | std::string
    * interval        | int
    * turn off        | bool
    * interpolate     | Teuchos::Array<std::string>
    *
    * parameter descriptions:
    * - out file, the base name to dump the visualization output.
    * - interval, dump the ouptut only at every \f$ i^{th} \f$ specified
    * interval for which \ref Output::write is called.
    * - turn off, completely turn off dumping output to file.
    * - interpolate, if this parameter is specified each field
    * corresponding to the names in this array will be interpolated onto
    * the mesh for output.
    *
    * \param d The relevant \ref goal::Discretization. */
  Output(RCP<const ParameterList> p, RCP<Discretization> d);

  /** \brief Write the visualization output.
    * \param t The time associated with the current output. */
  void write(const double t);

 private:
  void write_vtk(const double t);
  RCP<const ParameterList> params;
  RCP<Discretization> disc;
  Teuchos::Array<std::string> fields;
  bool turn_off;
  int interval;
  std::string name;
  int pos;
  int index;
};

}  // namesapce goal

#endif
