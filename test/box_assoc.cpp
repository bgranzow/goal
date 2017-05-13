#include <fstream>

namespace test {

void write_2D_assoc() {
  std::fstream assoc;
  assoc.open("box2D.txt", std::ios::out);
  assoc << "element set box 1\n"
    << "2 0\n"
    << "node set xmin 1\n"
    << "1 1\n"
    << "node set ymin 1\n"
    << "1 0\n"
    << "node set xmax 1\n"
    << "1 2\n"
    << "node set ymax 1\n"
    << "1 3\n"
    << "side set xmin 1\n"
    << "1 1\n"
    << "side set ymin 1\n"
    << "1 0\n"
    << "side set xmax 1\n"
    << "1 2\n"
    << "side set ymax 1\n"
    << "1 3\n";
  assoc.close();
}

void write_3D_assoc() {
  std::fstream assoc;
  assoc.open("box3D.txt", std::ios::out);
  assoc << "element set box 1\n"
    << "3 0\n"
    << "node set xmin 1\n"
    << "2 2\n"
    << "node set ymin 1\n"
    << "2 1\n"
    << "node set zmin 1\n"
    << "2 0\n"
    << "node set xmax 1\n"
    << "2 3\n"
    << "node set ymax 1\n"
    << "2 4\n"
    << "node set zmax 1\n"
    << "2 5\n"
    << "side set xmin 1\n"
    << "2 2\n"
    << "side set ymin 1\n"
    << "2 1\n"
    << "side set zmin 1\n"
    << "2 0\n"
    << "side set xmax 1\n"
    << "2 3\n"
    << "side set ymax 1\n"
    << "2 4\n"
    << "side set zmax 1\n"
    << "2 5\n";

  assoc.close();
}

} // end namespace test

int main() {
  test::write_2D_assoc();
  test::write_3D_assoc();
}
