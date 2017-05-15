#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <Teuchos_ParameterList.hpp>

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: disc");
  GOAL_ALWAYS_ASSERT(argc == 4);
  Teuchos::ParameterList p;
  p.set<std::string>("geom file", argv[1]);
  p.set<std::string>("mesh file", argv[2]);
  p.set<std::string>("assoc file", argv[3]);
  p.set<bool>("reorder mesh", true);
  p.set<int>("workset size", 25);
  auto d = goal::create_disc(p);
  goal::print("num dims: %d", d->get_num_dims());
  goal::print("ws size: %d", d->get_ws_size());
  goal::print("num elem sets: %d", d->get_num_elem_sets());
  goal::print("num side sets: %d", d->get_num_side_sets());
  goal::print("num node sets: %d", d->get_num_node_sets());
  for (int i = 0; i < d->get_num_elem_sets(); ++i) {
    goal::print("elem set: %d", i);
    goal::print(" name: %s", (d->get_elem_set_name(i)).c_str());
    goal::print(" elem type: %d", d->get_elem_type(i));
    goal::print(" num elem worksets: %d", d->get_num_elem_worksets(i));
  }
  for (int i = 0; i < d->get_num_side_sets(); ++i) {
    goal::print("side set: %d", i);
    goal::print(" name: %s", (d->get_side_set_name(i)).c_str());
    goal::print(" side type: %d", d->get_side_type(i));
    goal::print(" num side worksets: %d", d->get_num_side_worksets(i));
  }
  for (int i = 0; i < d->get_num_elem_sets(); ++i)
    for (int ws = 0; ws < d->get_num_elem_worksets(i); ++ws)
      d->get_elems(d->get_elem_set_name(i), ws);
  for (int i = 0; i < d->get_num_side_sets(); ++i)
    for (int ws = 0; ws < d->get_num_side_worksets(i); ++ws)
      d->get_sides(d->get_side_set_name(i), ws);
  d->update();
  goal::destroy_disc(d);
  goal::finalize();
}
