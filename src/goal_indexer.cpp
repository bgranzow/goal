#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

#include "goal_control.hpp"
#include "goal_discretization.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"

namespace goal {

using Teuchos::rcp;

Indexer::Indexer(RCP<Discretization> d, std::vector<RCP<Field> > f) {
  disc = d;
  fields = f;
  assert(fields.size() > 0);
  for (size_t i = 0; i < fields.size(); ++i)
    assert(fields[i]->get_associated_dof_idx() == int(i));
  comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  set_elem_block(0);
}

Indexer::~Indexer() {
  owned_map = Teuchos::null;
  ghost_map = Teuchos::null;
  owned_graph = Teuchos::null;
  ghost_graph = Teuchos::null;
  node_map = Teuchos::null;
  coords = Teuchos::null;
}

void Indexer::set_elem_block(const int idx) {
  elem_block_idx = idx;
  for (int i = 0; i < get_num_dof_fields(); ++i) fields[i]->set_elem_block(idx);
}

int Indexer::get_num_dof_fields() const { return fields.size(); }

int Indexer::get_num_total_elem_dofs() const {
  int dofs = 0;
  for (int f = 0; f < get_num_dof_fields(); ++f)
    dofs += fields[f]->get_num_elem_dofs();
  return dofs;
}

int Indexer::get_dof_idx(std::string const& name) const {
  int idx = -1;
  for (int f = 0; f < get_num_dof_fields(); ++f) {
    auto n = apf::getName(fields[f]->get_apf_field());
    if (n == name) idx = f;
  }
  if (idx < 0) fail("could not find dof field %s", name.c_str());
  return idx;
}

int Indexer::get_num_node_sets() const { return disc->get_num_node_sets(); }

std::vector<apf::Node> const& Indexer::get_node_set_nodes(
    std::string const& set, const int idx) {
  if (!node_sets.count(set)) fail("node set %s not found", set.c_str());
  return node_sets[set][idx];
}

apf::Mesh* Indexer::get_apf_mesh() { return disc->get_apf_mesh(); }

}  /* namespace goal */
