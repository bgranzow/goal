#ifndef goal_states_hpp
#define goal_states_hpp

/// @file goal_states.hpp

#include <MiniTensor.h>

/// @cond
namespace apf {
class Mesh;
class Field;
class MeshEntity;
}
/// @endcond

namespace goal {

using Teuchos::RCP;

/// @cond
class Discretization;
/// @endcond

/// @brief An interface for storing and querying state data.
/// @details This structure attempts to make it easier to manage and query
/// statie information defined at integration points.
class States {

  public:

    /// @brief Construct the state object.
    /// @param d The relevant \ref goal::Discretization.
    /// @param q The quadrature degree to use.
    States(Discretization* d, int q);

    /// @brief Destroy the state object.
    /// @details This will delete all APF fields built.
    ~States();

    /// @brief Add a state variable to this container.
    /// @param n The name of the state field.
    /// @param t The APF ValueType of hte field.
    /// @param s Should the state at the previous time be saved?
    /// @param I Should the state variable be initialized to the identity?
    void add(const char* n, int t, bool s = false, bool I = false);

    /// @brief Project the state fields to a different quadrature degree.
    /// @param q The degree to project the state to.
    void project(int q);

    /// @brief Update the state variables.
    /// @details Set values of old states to the current values.
    void update();

    /// @brief Set the value of a scalar at an integration point.
    /// @param n The name of the state variable.
    /// @param e The mesh entity the integration point is associated with.
    /// @param ip The integration point index in the mesh entity.
    /// @param v The value of the scalar to be set.
    template <typename T>
      void set_scalar(const char* n, apf::MeshEntity* e, int ip, T const& v);

    /// @brief Set the value of a vector at an integration point.
    /// @param n The name of the state variable.
    /// @param e The mesh entity the integration point is associated with.
    /// @param ip The integration point index in the mesh entity.
    /// @param v The value of the vector to be set. */
    template <typename T>
      void set_vector(const char* n, apf::MeshEntity* e, int ip,
          minitensor::Vector<T> const& v);

    /// \brief Set the value of a tensor at an integration point.
    /// \param n The name of the state variable.
    /// \param e The mesh entity the integration point is associated with.
    /// \param ip The integration point index in the mesh entity.
    /// \param v The value of the tensor to be set.
    template <typename T>
      void set_tensor(const char* n, apf::MeshEntity* e, int ip,
          minitensor::Tensor<T> const& v);

    /// \brief Get the value of a scalar at an integration point.
    /// \param n The name of the state variable.
    /// \param e The mesh entity the integration point is associated with.
    /// \param ip The integration point index in the mesh entity.
    /// \param v The return value of the scalar at the integration point.
    template <typename T>
      void get_scalar(const char* n, apf::MeshEntity* e, int ip, T& v);

    /// @brief Get the value of a vector at an integration point.
    /// @param n The name of the state variable.
    /// @param e The mesh entity the integration point is associated with.
    /// @param ip The integration point index in the mesh entity.
    /// @param v The return value of the vector at the integration point.
    template <typename T>
      void get_vector(const char* n, apf::MeshEntity* e, int ip,
          minitensor::Vector<T>& v);

    /// @brief Get the value of a tensor at an integration point.
    /// @param n The name of the state variable.
    /// @param e The mesh entity the integration point is associated with.
    /// @param ip The integration point index in the mesh entity.
    /// @param v The return value of the tensor at the integration point.
    template <typename T>
      void get_tensor(const char* n, apf::MeshEntity* e, int ip,
          minitensor::Tensor<T>& v);

  private:

    void add_state(const char* n, int t);
    void add_old_state(const char* n, int t);

    int q_degree;
    Discretization* disc;
    apf::Mesh* apf_mesh;
    std::vector<apf::Field*> states;
    std::vector<apf::Field*> old_states;

};

/// @brief Create a \ref goal::States object.
/// @param d The relevant \ref goal::Discretization.
/// @param q The quadrature degree to use.
States* create_states(Discretization* d, int q);

/// @brief Destroy a \ref goal::States object.
/// @param s The states object to destroy.
void destroy_states(States* s);

} // end namespace goal

#endif
