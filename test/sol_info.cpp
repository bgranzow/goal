#include <goal_control.hpp>
#include <goal_disc.hpp>
#include <goal_sol_info.hpp>

namespace test {

static void check_ops(goal::SolInfo* s) {
  s->resume_fill();
  s->zero_all();
  s->gather_all();
  s->complete_fill();
}

static void check_owned(goal::SolInfo* s) {
  auto R = s->owned->R;
  auto dMdu = s->owned->dMdu;
  auto dRdu = s->owned->dRdu;
  goal::print("owned R rows: %lu", R->getGlobalLength());
  goal::print("owned dMdu rows: %lu", dMdu->getGlobalLength());
  goal::print("owned dRdu rows: %lu", dRdu->getGlobalNumRows());
  goal::print("owned dRdu cols: %lu", dRdu->getGlobalNumCols());
}

static void check_ghost(goal::SolInfo* s) {
  auto R = s->ghost->R;
  auto dMdu = s->ghost->dMdu;
  auto dRdu = s->ghost->dRdu;
  goal::print("ghost R rows: %lu", R->getGlobalLength());
  goal::print("ghost dMdu rows: %lu", dMdu->getGlobalLength());
  goal::print("ghost dRdu rows: %lu", dRdu->getGlobalNumRows());
  goal::print("ghost dRdu cols: %lu", dRdu->getGlobalNumCols());
}

static void check_sol_info(goal::SolInfo* s) {
  check_ops(s);
  check_owned(s);
  check_ghost(s);
}

}

int main(int argc, char** argv) {
  goal::initialize();
  goal::print("unit test: sol info");
  GOAL_ALWAYS_ASSERT(argc == 4);
  Teuchos::ParameterList p;
  p.set<std::string>("geom file", argv[1]);
  p.set<std::string>("mesh file", argv[2]);
  p.set<std::string>("assoc file", argv[3]);
  auto d = goal::create_disc(p);
  d->build_data();
  auto s = goal::create_sol_info(d);
  test::check_sol_info(s);
  goal::destroy_sol_info(s);
  goal::destroy_disc(d);
  goal::finalize();
}
