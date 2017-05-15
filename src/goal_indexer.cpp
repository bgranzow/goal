#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"
#include "goal_strided_indexer.hpp"

namespace goal {

using Teuchos::rcp;

Indexer::Indexer(Discretization* d, std::vector<Field*> const& f) {
  disc = d;
  fields = f;
  GOAL_ALWAYS_ASSERT(fields.size() > 0);
  for (size_t i = 0; i < fields.size(); ++i)
    GOAL_ALWAYS_ASSERT(fields[i]->get_associated_dof_idx() == (int)i);
  comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
}

Indexer::~Indexer() {
  comm = Teuchos::null;
  owned_map = Teuchos::null;
  ghost_map = Teuchos::null;
  node_map = Teuchos::null;
  coords = Teuchos::null;
  owned_graph = Teuchos::null;
  ghost_graph = Teuchos::null;
}

apf::Mesh* Indexer::get_apf_mesh() {
  return disc->get_apf_mesh();
}

int Indexer::get_num_fields() const {
  return fields.size();
}

int Indexer::get_num_total_dofs(const int t) const {
  int dofs = 0;
  for (int f = 0; f < get_num_fields(); ++f)
    dofs += fields[f]->get_num_nodes(t);
  return dofs;
}

int Indexer::get_num_node_sets() const {
  return disc->get_num_node_sets();
}

std::vector<apf::Node> const& Indexer::get_node_set_nodes(
    std::string const& set, const int i) {
  if (! node_sets.count(set))
    fail("node set %s not found", set.c_str());
  return node_sets[set][i];
}

Indexer* create_indexer(
    int type, Discretization* d, std::vector<Field*> const& f) {
  Indexer* indexer = 0;
  if (type == STRIDED)
    indexer = new StridedIndexer(d, f);
  else
    fail("unknown indexer type");
  return indexer;
}

void destroy_indexer(Indexer* i) {
  delete i;
}

} // end namespace goal
