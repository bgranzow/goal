#ifndef goal_indexer_hpp
#define goal_indexer_hpp

/// @file goal_indexer.hpp

#include "goal_data_types.hpp"

/// @cond
namespace apf {
struct Node;
class Mesh;
class Field;
class FieldShape;
class MeshEntity;
template <class T> class NumberingOf;
typedef NumberingOf<int> Numbering;
typedef NumberingOf<long> GlobalNumbering;
}
/// @endcond

namespace goal {

/// @cond
class Field;
class Discretization;
/// @endcond

/// @brief The type of indexer.
enum IndexerType {
  /// @brief A strided indexer type.
  STRIDED = 0
};

/// @brief The DOF indexer interface.
/// @details The indexer structure is responsible for building the maps and
/// graphs that describe the parallel distribution of local to global data
/// in the linear algebra data structures. It does this by numbering the
/// DOFs associated with an input vector of DOF fields. This class is also
/// responsible for building node sets that are associated with the
/// specification of Dirichlet boundary conditions. This class is an
/// abstract base class, from which specific DOF indexing strategies
/// can be implemented.
///
/// See the Tpetra tutorial \cite TpetraTutorial for more details.
class Indexer {

  public:

    /// @brief Construct the indexer.
    /// @param d The relevant \ref goal::Discretization.
    /// @param f A vector of \ref goal::Field s to index.
    Indexer(Discretization* d, std::vector<Field*> const& f);

    /// @brief Destroy the indexer.
    /// @details This destroys maps, graphs, and multi-vectors that the
    /// indexer has constructed.
    virtual ~Indexer();

    /// @brief Get the owned DOF map.
    RCP<const Map> get_owned_map() { return owned_map; }

    /// @brief Get the ghost DOF map.
    RCP<const Map> get_ghost_map() { return ghost_map; }

    /// @brief Get the owned CRS element connectivity graph.
    RCP<const Graph> get_owned_graph() { return owned_graph; }

    /// @brief Get the ghost CRS element connectivity graph.
    RCP<const Graph> get_ghost_graph() { return ghost_graph; }

    /// @brief Get the coord multivector for owned nodes.
    /// @details Only non-null if all input fields have same basis type.
    RCP<MultiVector> get_coords() { return coords; }

    /// @brief Get the ith \ref goal::Field.
    Field* get_field(const int i) { return fields[i]; }

    /// @brief Get the underlying \ref goal::Discretization.
    Discretization* get_disc() { return disc; }

    /// @brief Get the underlying APF mesh.
    apf::Mesh* get_apf_mesh();

    /// @brief Get the total number of fields.
    int get_num_fields() const;

    /// @brief Get the total number of DOFs on this entity type.
    int get_num_total_dofs(const int t) const;

    /// @brief Get the number of node sets.
    int get_num_node_sets() const;

    /// @brief Get the nodes in a given node set.
    /// @param set The relevant node set name.
    /// @param i The relevant \ref goal::Field index.
    std::vector<apf::Node> const& get_node_set_nodes(
        std::string const& set, const int i);

    /// @brief Get the offset into elemental DOFs.
    /// @param i The \ref goal::Field index.
    /// @param n The elemental node associated with the field.
    virtual int get_elem_dof_offset(const int i, const int n) = 0;

    /// @brief Get the ghost local row ID.
    /// @param i The \ref goal::Field index.
    /// @param e The relevant mesh entity.
    /// @param n The local node associated with the mesh entity.
    virtual LO get_ghost_lid(
        const int i, apf::MeshEntity* e, const int n) = 0;

    /// @brief Get the owned local row ID.
    /// @param i The \ref goal::Field index.
    /// @param n An APF node data structure.
    virtual LO get_owned_lid(const int i, apf::Node const& n) = 0;

    /// @brief Get entity ghost local row IDs.
    /// @param e The relevant mesh entity.
    /// @param lids The returned local IDs.
    virtual void get_ghost_lids(
        apf::MeshEntity* e, std::vector<LO>& lids) = 0;

    /// @brief Sum values from the vector du into the fields f.
    /// @param f The \ref goal::Field s to sum into.
    /// @param du The relevant increment vector.
    virtual void add_to_fields(
        std::vector<Field*> const& f, RCP<Vector> du) = 0;

    /// @brief Set values from the vector x to the fields f.
    /// @param f The \ref goal::Field s to fill in.
    /// @param x The relevant vector to grab data from.
    virtual void set_to_fields(
        std::vector<Field*> const& f, RCP<Vector> x) = 0;

  protected:

    Discretization* disc;
    std::vector<Field*> fields;
    std::map<std::string,
      std::vector<std::vector<apf::Node> > > node_sets;

    RCP<const Comm> comm;
    RCP<const Map> owned_map;
    RCP<const Map> ghost_map;
    RCP<const Map> node_map;
    RCP<MultiVector> coords;
    RCP<Graph> owned_graph;
    RCP<Graph> ghost_graph;
};

/// @brief Create an indexer object.
/// @param type The \ref goal::IndexerType.
/// @param d The relevant \ref goal::Discretization.
/// @param f The \ref goal::Field s to index.
Indexer* create_indexer(
    int type, Discretization* d, std::vector<Field*> const& f);

/// @brief Destroy an indexer.
/// @param d The \ref goal::Indexer to destroy.
void destroy_indexer(Indexer* i);

} // end namespace goal

#endif
