#ifndef GOAL_STRIDED_INDEXER_HPP
#define GOAL_STRIDED_INDEXER_HPP

/** \file goal_strided_indexer.hpp */

#include "goal_indexer.hpp"

namespace goal {

/** \brief A strided DOF indexer class.
  * \details This class constructs a global indexing scheme using a strided
  * numbering of the degree of freedoms. For this indexer, all DOF \ref
  * goal::Field s must have the same basis type and be of the same
  * polynomial order. Each node in the mesh is numbered with this indexer
  * and then the DOF numbering is computed as:
  * - DOF = node_id * num_eqs + eq
  * where num_eqs is the total number of PDE equations (as determined by
  * the input DOF \ref goal::Field s) and eq is the current PDE equation
  * of interest. */
class StridedIndexer : public Indexer {
 public:
  /** \brief The strided indexer constructor.
    * \param d The relevant \ref goal::Discretization structure.
    * \param f A vector of DOF \ref goal::Field s. */
  StridedIndexer(RCP<Discretization> d, std::vector<RCP<Field> > f);

  /** \brief Destroy the strided indexer.
    * \details This will destroy all APF numberings associated with this
    * indexer. */
  ~StridedIndexer();

  /** \brief Get the offset into the elemental DOFs.
    * \param idx The \ref goal::Field DOF index.
    * \param n The elemental node associated with the field.
    * \param c The field component associated with the node. */
  int get_elem_dof_offset(const int idx, const int n, const int c);

  /** \brief Get the ghost local row ID.
    * \param idx The \ref goal::Field DOF index.
    * \param e The relevant mesh entity.
    * \param n The local node associated with the mesh entity.
    * \param c the field component associated with the local node. */
  LO get_ghost_lid(
      const int idx, apf::MeshEntity* e, const int n, const int c);

  /** \brief Get the owned local row ID.
    * \param idx The ref goal::Field DOF index.
    * \param n An APF node data structure (MeshEntity* + local node).
    * \param c The field component associated with the APF node. */
  LO get_owned_lid(const int idx, apf::Node const& n, const int c);

  /** \brief Get the elemental ghost local row IDs.
    * \param e The relevant mesh entity.
    * \param lids The returned local IDs. */
  void get_ghost_lids(apf::MeshEntity* e, std::vector<LO>& lids);

  /** \brief Sum values from the vector du into the fields f.
    * \param f The \ref goal::Field s to sum into.
    * \param du The relevant increment vector. */
  void add_to_fields(std::vector<RCP<Field> > const& f, RCP<Vector> du);

  /** \brief Set values from the vector x to the fields f.
    * \param f The \ref goal::Field s to fill in.
    * \param x The relevant solution vector. */
  void set_to_fields(std::vector<RCP<Field> > const& f, RCP<Vector> x);

 private:
  int num_eqs;
  void compute_owned_maps();
  void compute_ghost_map();
  void compute_graphs();
  void compute_coords();
  void compute_node_sets();
  apf::Numbering* owned_numbering;
  apf::Numbering* ghost_numbering;
  apf::GlobalNumbering* global_numbering;
};

} /* namespace goal */

#endif
