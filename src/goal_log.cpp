#include <sstream>
#include <cassert>
#include <iomanip>

#include "goal_control.hpp"
#include "goal_log.hpp"

namespace goal {

Log::Log(bool dual, bool exact, double J_exact) {
  J = J_exact;
  have_dual = dual;
  have_exact_J = exact;
}

void Log::print_banner() {
  std::ostringstream oss;
  oss << "time  iter  pDOFs  Ju_h";
  if (have_exact_J)
    oss << "  J   E";
  if (have_dual)
    oss << "  dDOFs  E_h  B_h";
  if (have_dual && have_exact_J)
    oss << "  I  IB";
  std::string banner = oss.str();
  print("%s", banner.c_str());
}

void Log::pre_validate() {
  std::size_t n = time.size();
  assert(iter.size() == n);
  assert(pDOFs.size() == n);
  assert(Ju_h.size() == n);
  if (have_dual) {
    assert(dDOFs.size() == n);
    assert(E_h.size() == n);
    assert(B_h.size() == n);
  }
}

void Log::compute_error(const int i) {
  if (! have_exact_J) return;
  E = J - Ju_h[i];
}

void Log::compute_effectivities(const int i) {
  if ((! have_exact_J) && (! have_dual)) return;
  I = E_h[i] / E;
  IB = B_h[i] / E;
}

void Log::print_index(const int i) {
  pre_validate();
  compute_error(i);
  compute_effectivities(i);
  std::ostringstream oss;
  oss << std::setprecision(10);
  oss << time[i] << "  ";
  oss << iter[i] << "  ";
  oss << pDOFs[i] << "  ";
  oss << Ju_h[i] << "  ";
  if (have_exact_J) {
    oss << J << "  ";
    oss << E << "  ";
  }
  if (have_dual) {
    oss << dDOFs[i] << "  ";
    oss << E_h[i] << "  ";
    oss << B_h[i] << "  ";
  }
  if (have_dual && have_exact_J) {
    oss << I << "  ";
    oss << IB << "  ";
  }
  std::string line = oss.str();
  print("%s", line.c_str());
}

void Log::print_current() {
  auto i = time.size()-1;
  print_banner();
  print_index(i);
}

void Log::print_summary() {
  print_banner();
  for (size_t i = 0; i < time.size(); ++i)
    print_index(i);
}

} /* namespace goal */
