#include <PCU.h>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "goal_control.hpp"
#include "goal_log.hpp"

namespace goal {

Log::Log(bool dual, bool exact, double J_exact) {
  J = J_exact;
  have_dual = dual;
  have_exact_J = exact;
}

std::string Log::get_banner() {
  std::ostringstream oss;
  oss << "time  iter  pDOFs  Ju_h";
  if (have_exact_J)
    oss << "  J   E";
  if (have_dual)
    oss << "  dDOFs  E_h  B_h";
  if (have_dual && have_exact_J)
    oss << "  I  IB";
  std::string banner = oss.str();
  return banner;
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

std::string Log::get_index(const int i) {
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
    oss << std::abs(E) << "  ";
  }
  if (have_dual) {
    oss << dDOFs[i] << "  ";
    oss << std::abs(E_h[i]) << "  ";
    oss << std::abs(B_h[i]) << "  ";
  }
  if (have_dual && have_exact_J) {
    oss << std::abs(I) << "  ";
    oss << std::abs(IB) << "  ";
  }
  std::string line = oss.str();
  return line;
}

void Log::print_index(const int i) {
  pre_validate();
  print("%s", get_index(i).c_str());
}

void Log::print_current() {
  auto i = time.size()-1;
  print("%s", get_banner().c_str());
  print_index(i);
}

void Log::print_summary() {
  print("%s", get_banner().c_str());
  for (size_t i = 0; i < time.size(); ++i)
    print_index(i);
}

void Log::print_summary(std::string const& n) {
  if (PCU_Comm_Self()) return;
  std::fstream fs;
  fs.open(n.c_str(), std::ios::out);
  fs << get_banner() << std::endl;
  for (size_t i = 0; i < time.size(); ++i)
    fs << get_index(i) << std::endl;
  fs.close();
}

} /* namespace goal */
