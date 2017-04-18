#ifndef GOAL_FIELD_HPP
#define GOAL_FIELD_HPP

/** \file goal_field.hpp */

#include <Teuchos_RCP.hpp>

#include "goal_field_types.hpp"

/** \cond */
namespace PHX {
class DataLayout;
}

namespace apf {
class Mesh;
class Field;
class FieldShape;
}
/** \endcond */

namespace goal {

using Teuchos::RCP;

/** \cond */
class Discretization;
/** \endcond */

/** \brief Field information structure.
  * \details This structure contains all of the necessary information to
  * fully define a \ref goal::Field. */
struct FieldInfo {
  /** \brief The \ref goal::Discretization structure */
  RCP<Discretization> disc;

  /** \brief The unique name of the field. */
  std::string name;

  /** \brief The polynomial order of the field. */
  int p_order;

  /** \brief The quadrature degree used for the field. */
  int q_degree;

  /** \brief The \ref goal::FieldType of this field. */
  int value_type;

  /** \brief The \ref goal::BasisType of this field. */
  int basis_type;
};

/** \brief A container for field information.
  * \details The \ref goal::Field container augments an APF field with
  * naming conventions for outputs that arise from operations on this
  * field, as well element-level dimension information that is specific
  * to the field. */
class Field {
 public:
  /** \brief Construct a field.
    * \details The \ref goal::FieldInfo struct that defines the field. */
  Field(FieldInfo* info);

  /** \brief Destroy the field.
    * \details This will destroy the underlying APF field. */
  ~Field();

  /** \brief Set the current element block for this field.
    * \details This is useful for scenarios where there are different
    * element topologies in different element blocks, so that the
    * element-level dimension information provided by the \ref goal::Field
    * APIs is correct. */
  void set_elem_block(const int idx);

  /** \brief Set the derivative seed of this field.
    * \details If the field is a degree of freedom, this should be 1. If the
    * field is a linear combination of the degrees of freedom (e.g. for
    * time-derivative fields computed by finite-differencing) it should be
    * set to an appropriate scaling factor (e.g. 1/dt). If the field is not
    * a degree of freedom, it should be set to 0. */
  void set_seed(double s) { seed = s; }

  /** \brief Set the degree of freedom status of this field. */
  void set_dof_status(bool is) { dof = is; }

  /** \brief Associate this field with a particular DOF index.
    * \details The associated DOF index of a field is defaulted to -1,
    * indicating that the field is not associated with any particular DOF
    * variable. Fields that are DOF fields should be associated with an
    * appropriate index. Other fields, such as time derivative fields, that
    * are not actual DOF fields, should still associated with a DOF index. */
  void set_associated_dof_idx(const int i) { idx = i; }

  /** \brief Is this field a degree of freedom field? */
  bool is_dof() { return dof; }

  /** \brief Get the underlying APF mesh data structure.
    * \details This API is mainly provided for convenience. */
  apf::Mesh* get_apf_mesh() const { return apf_mesh; }

  /** \brief Get the underlying APF field data structure. */
  apf::Field* get_apf_field() const { return apf_field; }

  /** \brief Get the underlying APF field shape.
    * \details This data fully defines the basis functions used for this
    * field. */
  apf::FieldShape* get_apf_basis() const { return apf_basis; }

  /** \brief Returns the seed value for this field. */
  double get_seed_value() const { return seed; }

  /** \brief Get the polynomial order of this field's basis functions. */
  int get_p_order() const { return p_order; }

  /** \brief Get the quadrature degree used for this field. */
  int get_q_degree() const { return q_degree; }

  /** \brief Get the associated dof index for this field.
    * \details See \ref goal::Field::set_associated_dof_idx. */
  int get_associated_dof_idx() const { return idx; }

  /** \brief Get the \ref goal::FieldType of this field. */
  int get_value_type() const { return value_type; }

  /** \brief Get the \ref goal::BasisType of this field. */
  int get_basis_type() const { return basis_type; }

  /** \brief Get the number of spatial dimensions for this field. */
  int get_num_dims() const;

  /** \brief Get this field's number of components.
    * \details
    * - \ref goal::SCALAR = 1
    * - \ref goal::VECTOR = \ref goal::Field::get_num_dims. */
  int get_num_components() const;

  /** \brief Get this field's number of elemental nodes. */
  int get_num_elem_nodes() const;

  /** \brief Get this field's number of elemental vertices. */
  int get_num_elem_vtx() const;

  /** \brief Get this field's number of elemental DOFs.
    * \details This is equal to the number of elemental nodes times the
    * number of field components. */
  int get_num_elem_dofs() const;

  /** \brief Get this field's number of elemental integration points. */
  int get_num_elem_ips() const;

  /** \brief Returns this field's name. */
  std::string get_name() const;

  /** \brief Returns this field's gradient name.
    * \details = "grad_" + \ref goal::Field::get_name. */
  std::string get_grad_name() const;

  /** \brief Returns this field's basis name.
    * \details = \ref goal::BasisType + \ref goal::Field::get_p_order.
    *
    * (e.g. LAGRANGE2) */
  std::string get_basis_name() const;

  /** \brief Returns this field's basis gradient name.
    * \details = "grad_" + \ref goal::Field::get_basis_name. */
  std::string get_grad_basis_name() const;

  /** \brief Get this field's residual name.
    * \details = \ref goal::Field::get_name + "_residual". */
  std::string get_residual_name() const;

  /** \brief Get this field's weighted differential volume name.
    * \details = \ref goal::Field::get_name + "_wdv". */
  std::string get_wdv_name() const;

  /** \brief Get the data layout.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, Node)
    * \ref goal::VECTOR     | (Elem, Node, Dim) */
  RCP<PHX::DataLayout> get_dl();

  /** \brief Get the gradient data layout.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, Node, Dim)
    * \ref goal::VECTOR     | (Elem, Node, Dim, Dim) */
  RCP<PHX::DataLayout> get_grad_dl();

  /** \brief Get the basis weight data layout.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, Node, IP)
    * \ref goal::VECTOR     | (Elem, Node, IP, Dim) */
  RCP<PHX::DataLayout> get_weight_dl();

  /** \brief Get the basis gradient weight data layout.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, Node, IP, Dim)
    * \ref goal::VECTOR     | (Elem, Node, IP, Dim, Dim) */
  RCP<PHX::DataLayout> get_grad_weight_dl();

  /** \brief Get the field interpolated to IPs data layout.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, IP)
    * \ref goal::VECTOR     | (Elem, IP, Dim) */
  RCP<PHX::DataLayout> get_interpolated_dl();

  /** \brief Get the gradient of the field interpolated to IPs data layout.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, IP, Dim)
    * \ref goal::VECTOR     | (Elem, IP, Dim, Dim) */
  RCP<PHX::DataLayout> get_grad_interpolated_dl();

  /** \brief Get the data layout of a partition of unity.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, Node, IP)
    * \ref goal::VECTOR     | (Elem, Node, IP, Dim) */
  RCP<PHX::DataLayout> get_PU_dl();

  /** \brief Get the data layout of the gradient of a partition of unity.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, Node, IP, Dim)
    * \ref goal::VECTOR     | (Elem, Node, IP, Dim, Dim) */
  RCP<PHX::DataLayout> get_grad_PU_dl();

  /** \brief Get the data layout of the residual of a partition of unity.
    * \details
    * \ref goal::FieldType  | data layout
    * --------------------  | -----------
    * \ref goal::SCALAR     | (Elem, Node)
    * \ref goal::VECTOR     | (Elem, Node, Dim) */
  RCP<PHX::DataLayout> get_residual_PU_dl();

  /** \brief Get the data layout of a scalar at IPs.
    * \details data layout: (Elem, IP) */
  RCP<PHX::DataLayout> get_scalar_ip_dl();

  /** \brief Get the data layout of a vector at IPs.
    * \details data layout: (Elem, IP, Dim) */
  RCP<PHX::DataLayout> get_vector_ip_dl();

  /** \brief Get the data layout of a tensor at IPs.
    * \details data layout: (Elem, IP, Dim, Dim) */
  RCP<PHX::DataLayout> get_tensor_ip_dl();

  /** \brief Get the data layout of a scalar associated with elements.
    * \details data layout: (Elem) */
  RCP<PHX::DataLayout> get_elem_scalar_dl();

 private:
  RCP<Discretization> disc;
  apf::Mesh* apf_mesh;
  apf::Field* apf_field;
  apf::FieldShape* apf_basis;
  std::string name;
  int p_order;
  int q_degree;
  int value_type;
  int basis_type;
  int idx;
  bool dof;
  double seed;
  int elem_type;
};

/** \brief Project data from one field to another. */
void project_field(RCP<Field> to, RCP<Field> from);

}  /* namespace goal */

#endif
