#ifndef GOAL_STATE_FIELDS_HPP
#define GOAL_STATE_FIELDS_HPP

/** \file goal_state_fields.hpp */

#include <MiniTensor.h>

#include "goal_field_types.hpp"

/** \cond */
namespace apf {
class Mesh;
class Field;
class MeshEntity;
}
/** \endcond */

namespace goal {

using Teuchos::RCP;

/** \cond */
class Discretization;
/** \endcond */

/** \brief The state fields interface.
  * \details This structure attempts to make it easier to manage storing
  * and querying of state information defined at integration points in
  * applications that are developed using goal. */
class StateFields {
 public:
  /** \brief Construct a state fields structure.
    * \param d The relevant \ref goal::Discretization.
    * \param q_degree The quadrature degree to use for these states. */
  StateFields(RCP<Discretization> d, int q_degree);

  /** \brief Add a state variable to this container.
    * \param n The name of the state field.
    * \param type The \ref goal::FieldType of the field.
    * \param save Should the state at the previous time be saved?
    * \param I Should the state variable be initialized to the identity? */
  void add(const char* n, int type, bool save = false, bool I = false);

  /** \brief Project all state fields to a different quadrature degree. */
  void project(int q_degree);

  /** \brief Update the state variables.
    * \details Set the value of all states initialized with save=true at
    * \f$ t_{n-1} \f$ to be the newly computed state at \f$ t_n \f$. */
  void update();

  /** \brief Set the value of a scalar at an integration point.
    * \param name The name of the state variable.
    * \param e The mesh entity the integration point is associated with.
    * \param ip The integration point index in the mesh entity.
    * \param v The value of the scalar to be set. */
  template <typename T>
  void set_scalar(const char* name, apf::MeshEntity* e, int ip, T const& v);

  /** \brief Set the value of a vector at an integration point.
    * \param name The name of the state variable.
    * \param e The mesh entity the integration point is associated with.
    * \param ip The integration point index in the mesh entity.
    * \param v The value of the vector to be set. */
  template <typename T>
  void set_vector(const char* name, apf::MeshEntity* e, int ip,
      minitensor::Vector<T> const& v);

  /** \brief Set the value of a tensor at an integration point.
    * \param name The name of the state variable.
    * \param e The mesh entity the integration point is associated with.
    * \param ip The integration point index in the mesh entity.
    * \param v The value of the tensor to be set. */
  template <typename T>
  void set_tensor(const char* name, apf::MeshEntity* e, int ip,
      minitensor::Tensor<T> const& v);

  /** \brief Get the value of a scalar at an integration point.
    * \param name The name of the state variable.
    * \param e The mesh entity the integration point is associated with.
    * \param ip The integration point index in the mesh entity.
    * \param v The return value of the scalar at the integration point. */
  template <typename T>
  void get_scalar(const char* name, apf::MeshEntity* e, int ip, T& v);

  /** \brief Get the value of a vector at an integration point.
    * \param name The name of the state variable.
    * \param e The mesh entity the integration point is associated with.
    * \param ip The integration point index in the mesh entity.
    * \param v The return value of the vector at the integration point. */
  template <typename T>
  void get_vector(
      const char* name, apf::MeshEntity* e, int ip, minitensor::Vector<T>& v);

  /** \brief Get the value of a tensor at an integration point.
    * \param name The name of the state variable.
    * \param e The mesh entity the integration point is associated with.
    * \param ip The integration point index in the mesh entity.
    * \param v The return value of the tensor at the integration point. */
  template <typename T>
  void get_tensor(
      const char* name, apf::MeshEntity* e, int ip, minitensor::Tensor<T>& v);

 private:
  int q_degree;
  RCP<Discretization> disc;
  apf::Mesh* apf_mesh;
  std::vector<apf::Field*> states;
  std::vector<apf::Field*> old_states;
};

}  /* namespace goal */

#endif
