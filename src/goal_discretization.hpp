#ifndef GOAL_DISCRETIZATION_HPP
#define GOAL_DISCRETIZATION_HPP

/** \file goal_discretization.hpp */

#include <Teuchos_RCP.hpp>

/** \cond */
namespace Teuchos {
class ParameterList;
}

namespace apf {
struct Node;
struct StkModels;
class Mesh;
class Mesh2;
class MeshEntity;
}
/** \endcond */

namespace goal {

using Teuchos::RCP;
using Teuchos::ParameterList;

/** \brief The discretization interface.
  * \details The discretization object exists to augment an APF mesh data
  * structure with a geometric domain definition of the problem to be solved.
  * The discretization is taksed with computing element sets and facet sets,
  * which help define the domain definition of the problem. */
class Discretization {
 public:
  /** \brief The discretization constructor.
    * \param p A parameterlist fully describing the discretization.
    *
    * Parameter Name  | Parameter Type
    * --------------  | --------------
    * geom file       | std::string
    * mesh file       | std::string
    * assoc file      | std::string
    * reorder mesh    | bool
    * quadratic       | bool
    * workset size    | int
    * mesh            | apf::Mesh2*
    * associations    | apf::StkModels*
    *
    * \details If the parameter 'mesh' exists, then the Discretization
    * object will be constructed from an existing APF mesh, otherwise,
    * a mesh will be loaded from file using the 'geom file' and
    * 'mesh file' parameters. If the parameter 'associations' exists,
    * then the Discretization object will be constructed using an existing
    * associations data structure, otherwise, these associations will be
    * loaded from file using the 'assoc file' parameter. If the parameter
    * 'coord order' exists, then the discretization structure will attempt
    * to modify the APF mesh's underlying geometry representation to the
    * given order. The optional parameter 'quadratic' is used to modify
    * the geometry to a quadratic representation if it is not already. */
  Discretization(RCP<const ParameterList> p);

  /** \brief Discretization destructor.
    * \details This will destroy the APF mesh if it was loaded from file.
    * This will destroy the APF associtiations data if it was loaded from
    * file. */
  ~Discretization();

  /** \brief Update the discretization data structures.
    * \details This method is useful after the mesh topology has been
    * modified due to mesh adaptation. It calls in sequence:
    *
    * - compute_elem_blocks()
    * - compute_side_sets() */
  void update();

  /** \brief Returns the underlying APF mesh. */
  apf::Mesh2* get_apf_mesh() { return mesh; }

  /** \brief Returns the set definitions.
    * \details These sets are responsible for associating groups of domain
    * entities to provide a domain definition of the problem. */
  apf::StkModels* get_model_sets() { return sets; }

  /** \brief Returns the number of spatial dimensions of the mesh. */
  int get_num_dims() const { return num_dims; }

  /** \brief Returns the workset size.
    * \details This is the maximum number of mesh entities operated on
    * during a call to PHX::evaluateFields, and is used to limit memory
    * consumption and improve cache performance. */
  int get_ws_size() const { return ws_size; }

  /** \brief Return the number of element blocks. */
  int get_num_elem_blocks() const;

  /** \brief Returns the number of side sets. */
  int get_num_side_sets() const;

  /** \brief Returns the number of node sets. */
  int get_num_node_sets() const;

  /** \brief Returns the ith element block name. */
  std::string get_elem_block_name(const int idx) const;

  /** \brief Returns the ith side set name. */
  std::string get_side_set_name(const int idx) const;

  /** \brief Returns the ith node set name. */
  std::string get_node_set_name(const int idx) const;

  /** \brief Returns the element type associated with ith element block. */
  int get_elem_type(const int block_idx);

  /** \brief Returns the number of worksets in the ith element block. */
  int get_num_worksets(const int block_idx);

  /** \brief Returns the elements in a workset.
    * \param elem_block The element block to return elements from.
    * \param ws_idx The workset index in the relevant element block. */
  std::vector<apf::MeshEntity*> const& get_elems(
      std::string const& elem_block, const int ws_idx);

  /** \brief Return the facets associated with a side set.
    * \param side_set The name of the relevant side set. */
  std::vector<apf::MeshEntity*> const& get_sides(std::string const& side_set);

 private:
  void compute_elem_blocks();
  void compute_side_sets();
  RCP<const ParameterList> params;
  bool owns_mesh;
  bool owns_sets;
  int num_dims;
  int ws_size;
  apf::Mesh2* mesh;
  apf::StkModels* sets;
  std::map<std::string, std::vector<std::vector<apf::MeshEntity*> > >
      elem_blocks;
  std::map<std::string, std::vector<apf::MeshEntity*> > side_sets;
};

}  // namespace goal

#endif
