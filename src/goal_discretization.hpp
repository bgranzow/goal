#ifndef goal_discretization_hpp
#define goal_discretization_hpp

/// @file goal_discretization.hpp

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

/// @cond
namespace apf {
struct Node;
struct StkModels;
class Mesh2;
class MeshEntity;
}
/// @endcond

namespace goal {

using Teuchos::ParameterList;

/// @brief The discretization interface.
/// @details The discretization object exists to augment an APF mesh data
/// structure with a geometric domain definition of the problem to be solved.
/// This class is tasked with computing element sets and side sets, which
/// help define the domain definition of the problem.
class Discretization {

  public:

    /// @brief The discretization constructor.
    /// @param p A parameter list fully describing the discretization.
    /// parameter name  | parameter type
    /// --------------  | --------------
    /// geom file       | std::string
    /// mesh file       | std::string
    /// assoc file      | std::string
    /// reorder mesh    | bool
    /// make quadratic  | bool
    /// workset size    | int
    /// mesh            | apf::Mesh2*
    /// associations    | apf::StkModels*
    ///
    /// @details If the parameter 'mesh' exists, then the Discretization
    /// object will be constructed from an existing APF mesh. Otherwise,
    /// a mesh will be loaded from file using the 'geom file' and 'mesh file'
    /// parameters. If the parameter 'associations' exists, then the
    /// Discretization object will be constructed using an existing
    /// associations data structure. Otherwise, these associations will be
    /// loaded from file using the 'assoc file' parameter. If the parameter
    /// 'make quadratic' is set to true, then the discretization will attempt
    /// to modify the APF mesh's underlying geometry representation to a
    /// quadratic representation.
    Discretization(ParameterList const& p);

    /// @brief Discretization destructor.
    /// @details This will destroy the APF mesh if it was loaded from file
    /// and the APF associations data if it was loaded was from file.
    ~Discretization();

    /// @brief Update the discretization data structures.
    /// @details This method is useful after the mesh topology has been
    /// modiied due to mesh adaptation. It calls in sequence:
    /// - compute_elem_sets()
    /// - compute_side_sets()
    void update();

    /// @brief Returns the underlying APF mesh.
    apf::Mesh2* get_apf_mesh() { return mesh; }

    /// @brief Returns the set definitions.
    /// @details These sets are responsible for associating groups of
    /// domain entities to provide a domain definition of the problem.
    apf::StkModels* get_model_sets() { return sets; }

    /// @brief Returns the number of spatial dimensions of the mesh.
    int get_num_dims() const { return num_dims; }

    /// @brief Returns the workset size.
    /// @details This is the maximum number of mesh entities operated on
    /// during a call to PHX::evaluateFields, and is used to limit memory
    /// consumption and improve cache performance.
    int get_ws_size() const { return ws_size; }

    /// @brief Returns the number of element sets.
    int get_num_elem_sets() const;

    /// @brief Returns the number of side sets.
    int get_num_side_sets() const;

    /// @brief Returns the number of node sets.
    int get_num_node_sets() const;

    /// @brief Returns the ith element set name.
    std::string get_elem_set_name(const int i) const;

    /// @brief Returns the ith side set name.
    std::string get_side_set_name(const int i) const;

    /// @brief Returns the ith node set name.
    std::string get_node_set_name(const int i) const;

    /// @brief Returns the index of the element set of name n.
    int get_elem_set_idx(std::string const& n) const;

    /// @brief Returns the index of the side set of name n.
    int get_side_set_idx(std::string const& n) const;

    /// @brief Returns the index of the node set of name n.
    int get_node_set_idx(std::string const& n) const;

    /// @brief Returns the element entity type of the ith element set.
    int get_elem_type(const int i);

    /// @brief Returns the side entity type of the ith side set.
    int get_side_type(const int i);

    /// @brief Returns the number of worksets in the ith element set.
    int get_num_elem_worksets(const int i);

    /// @brief Returns the number of worksets in the ith side set.
    int get_num_side_worksets(const int i);

    /// @brief Returns the elements in a workset.
    /// @param elem_set The element set to return elements from.
    /// @param ws_idx The workset index in the relevant element set.
    std::vector<apf::MeshEntity*> const& get_elems(
        std::string const& elem_set, const int ws_idx);

    /// @brief Returns the sides in a workset.
    /// @param side_set The side set to return sides from.
    /// @param ws_idx The workset index in the relevant side set.
    std::vector<apf::MeshEntity*> const& get_sides(
        std::string const& side_set, const int ws_idx);

  private:

    void compute_elem_sets();
    void compute_side_sets();

    ParameterList params;

    bool owns_mesh;
    bool owns_sets;

    int num_dims;
    int ws_size;

    apf::Mesh2* mesh;
    apf::StkModels* sets;

    std::map<std::string,
      std::vector<std::vector<apf::MeshEntity*> > > elem_sets;
    std::map<std::string,
      std::vector<std::vector<apf::MeshEntity*> > > side_sets;
};


/// @brief Create a discretization.
/// @param p The discretization input parameters.
Discretization* create_disc(ParameterList const& p);

/// @brief Destroy a discretization.
/// @param d The discretization to destroy.
void destroy_disc(Discretization* d);

} // end namespace goal

#endif
