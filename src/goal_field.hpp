#ifndef goal_field_hpp
#define goal_field_hpp

/// @file goal_field.hpp

#include <Teuchos_RCP.hpp>

/// @cond
namespace PHX {
class DataLayout;
}

namespace apf {
class Mesh;
class Field;
class FieldShape;
}
/// @endcond

namespace goal {

using Teuchos::RCP;

/// @cond
class Discretization;
/// @endcond

/// @brief Available FEM basis types.
enum BasisType {
  /// @brief Lagrangian basis functions.
  LAGRANGE = 0,
  /// @brief Hierarchical basis functions.
  HIERARCHICAL = 1
};

/// @brief Field information structure
/// @details This structure contains all of the necessary
/// information to fully define a \ref goal::Field.
struct FieldInfo {
  /// @brief The relevant \ref goal::Discretization.
  Discretization* disc;
  /// @brief The unique name of the field.
  std::string name;
  /// @brief The polynomial order of the field.
  int p_order;
  /// @brief The quadrature degree used for the field.
  int q_degree;
  /// @brief The \ref goal::BasisType of this field.
  int basis_type;
};

/// @brief A container for field information.
/// @details The \ref goal::Field container augments an APF field with
/// naming conventions for outputs that arise form operations on this
/// field, as well as entity-level information that is specific to the field.
class Field {

  public:

    /// @brief Construct a field.
    /// @param i The \ref goal::FieldInfo that defines the field.
    Field(FieldInfo const& i);

    /// @brief Destroy the field.
    /// @details This will destroy the underlying APF field.
    ~Field();

    /// @brief Set the derivative seed of this field.
    /// @details Set the FAD seed value for this field.
    void set_seed(const double s) { seed = s; }

    /// @brief Associate this field with a DOF index.
    /// @details The DOF index is defaulted to -1.
    void set_associated_dof_idx(const int i) { idx = i; }

    /// @brief Returns the underlying APF mesh data structure.
    apf::Mesh* get_apf_mesh() const { return apf_mesh; }

    /// @brief Returns the underlying apf field structure.
    apf::Field* get_apf_field() const { return apf_field; }

    /// @brief Returns the underlying APF basis shape.
    apf::FieldShape* get_apf_basis() const { return apf_basis; }

    /// @brief Returns the polynomial order of this field.
    int get_p_order() const { return p_order; }

    /// @brief Returns the quadrature degree for this field.
    int get_q_degree() const { return q_degree; }

    /// @brief Get the \ref goal::BasisType of this field.
    int get_basis_type() { return basis_type; }

    /// @brief Returns the seed value for this field.
    double get_seed_value() { return seed; }

    /// @brief Returns the associated dof index for this field.
    int get_associated_dof_idx() { return idx; }

    /// @brief Returns the number of spatial dimensions.
    int get_num_dims() const;

    /// @brief Returns the number of entity nodes.
    /// @param t The entity type.
    int get_num_nodes(const int t) const;

    /// @brief Returns the number of entity vertices.
    /// @param t The entity type.
    int get_num_vtx(const int t) const;

    /// @brief Returns the number of entity integration points.
    /// @param t The entity type.
    int get_num_ips(const int t) const;

    /// @brief Returns the name of this field.
    std::string name() const;

    /// @brief Returns the name of the gradient of this field.
    /// @details "grad_" + name().
    std::string g_name() const;

    /// @brief Returns the name of the residual of this field.
    /// @details name() + "_resid".
    std::string resid_name() const;

    /// @brief Returns the name of the basis for this field.
    /// @details e.g. lagrange2
    std::string basis_name() const;

    /// @brief Returns the name of the gradient of the basis.
    /// @details "grad_" + basis_name().
    std::string g_basis_name() const;

    /// @brief Returns the name of the differential volume.
    /// @details basis_name + "_wdv".
    std::string wdv_name() const;

    /// @brief Returns the field data layout.
    /// @param t The entity type.
    /// @details (Ent, Node).
    RCP<PHX::DataLayout> dl(const int t);

    /// @brief Returns the gradient data layout.
    /// @details (Ent, Node, Dim).
    RCP<PHX::DataLayout> g_dl(const int t);

    /// @brief Returns the weight data layout.
    /// @details (Ent, Node, IP).
    RCP<PHX::DataLayout> w_dl(const int t);

    /// @brief Returns the gradient weight data layout.
    /// @details (Ent, Node, IP, Dim).
    RCP<PHX::DataLayout> g_w_dl(const int t);

    /// @brief Returns the data layout for field at ips.
    /// @details (Ent, IP)
    RCP<PHX::DataLayout> ip_dl(const int t);

    /// @brief Returns the data layout for grad field at ips.
    /// @details (Ent, IP, Dim).
    RCP<PHX::DataLayout> g_ip_dl(const int t);

    /// @brief Returns the data layout for a partition of unity.
    /// @details (Ent, Node, IP).
    RCP<PHX::DataLayout> PU_dl(const int t);

    /// @brief Returns the data layout for gradient of partiton of unity.
    /// @details (Ent, Node, IP, Dim).
    RCP<PHX::DataLayout> g_PU_dl(const int t);

    /// @brief Returns data layout for a scalar at ips.
    /// @details (Ent, IP).
    RCP<PHX::DataLayout> ip0_dl(const int t);

    /// @brief Returns data layout for a vector at ips.
    /// @details (Ent, IP, Dim).
    RCP<PHX::DataLayout> ip1_dl(const int t);

    /// @brief Returns data layout for a tensor at ips.
    /// @details (Ent, IP, Dim, Dim).
    RCP<PHX::DataLayout> ip2_dl(const int t);

    /// @brief Returns data layout for a scalar at the entity level.
    /// @details (Ent).
    RCP<PHX::DataLayout> ent0_dl(const int t);

  private:

    Discretization* disc;

    apf::Mesh* apf_mesh;
    apf::Field* apf_field;
    apf::FieldShape* apf_basis;

    std::string myname;
    int p_order;
    int q_degree;
    int basis_type;

    int idx;
    double seed;
};

/// @brief Create a field.
/// @param i The relevant \ref goal::FieldInfo data.
Field* create_field(FieldInfo const& i);

/// @brief Destroy a field.
/// @param f The field to destroy.
void destroy_field(Field* f);

} // end namespace goal

#endif
