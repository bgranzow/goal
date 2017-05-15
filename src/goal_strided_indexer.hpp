#ifndef goal_strided_indexer_hpp
#define goal_strided_indexer_hpp

/// @file goal_strided_indexer.hpp

#include "goal_indexer.hpp"

namespace goal {

/// @brief A strided DOF indexer class.
/// @details This class constructs a global indexing scheme using a
/// strided numbering of the degree of freedoms. For this indexer, all
/// \ref goal::Field s must have the same basis type and be of the same
/// polynomial order. Each node in the mesh is numbered with this indexer
/// and then the DOF numbering is computed as:
/// - DOF = node_id * num_eqs + eq
/// where num_eqs is the total number of PDE equations and eq is the
/// current PDE equation of interest.
class StridedIndexer : public Indexer {

  public:

    /// @brief Construct the strided indexer.
    /// @param d The relevant \ref goal::Discretization.
    /// @param f A vector of \ref goal::Field s to index.
    StridedIndexer(Discretization* d, std::vector<Field*> const& f);

    /// @brief Destroy the dstrided indexer.
    /// @details This destroys APF numberings built by this indexer.
    ~StridedIndexer();

    /// @brief Get the offset into the elemental DOFs.
    /// @param i The \ref goal::Field index.
    /// @param n The elemental node associated with the field.
    int get_elem_dof_offset(const int i, const int n);

    /// @brief Get the ghost local row ID.
    /// @param i The \ref goal::Field index.
    /// @param e The relevant mesh entity.
    /// @param n The local node associated with the mesh entity.
    LO get_ghost_lid(const int i, apf::MeshEntity* e, const int n);

    /// @brief Get the owned local row ID.
    /// @param i The \ref goal::Field index.
    /// @param n An APF node data structure.
    LO get_owned_lid(const int i, apf::Node const& n);

    /// @brief Get entity ghost local row IDs.
    /// @param e The relevant mesh entity.
    /// @param lids The returned local IDs.
    void get_ghost_lids(apf::MeshEntity* e, std::vector<LO>& lids);

    /// @brief Sum values from the vector du into the fields f.
    /// @param f The \ref goal::Field s to sum into.
    /// @param du The relevant increment vector.
    void add_to_fields(std::vector<Field*> const& f, RCP<Vector> du);

    /// @brief Set values from the vector x to the fields f.
    /// @param f The \ref goal::Field s to fill in.
    /// @param x The relevant vector to grab data from.
    void set_to_fields(std::vector<Field*> const& f, RCP<Vector> x);

  private:

    void compute_owned_maps();
    void compute_ghost_map();
    void compute_graphs();
    void compute_coords();
    void compute_node_sets();

    int num_eqs;
    apf::Mesh* mesh;
    apf::FieldShape* shape;
    apf::Numbering* owned_nmbr;
    apf::Numbering* ghost_nmbr;
    apf::GlobalNumbering* global_nmbr;
};

} // end namespace goal

#endif
