// Run the real Albert-Young case on the preprocessed diffusion data.
// Usage: albert_young_reference <converted.bin> <output.csv>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "albert_young_case.h"
#include "solver.h"

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " <converted.bin> <output.csv>\n";
    return 2;
  }
  albert_young::Case data(argv[1]);
  sayram2d::Solver solver(data.mesh());

  const std::size_t n1 = data.mesh().n1();
  const std::size_t n2 = data.mesh().n2();
  const double dt = 0.002;
  const int steps = 10;

  const std::vector<double> f0 = data.f();
  for (int step = 0; step < steps; ++step) {
    solver.solve(data.equation(), data.f_const_view(), data.f_view(), dt);
  }
  const std::vector<double>& f10 = data.f();

  std::ofstream out(argv[2]);
  if (!out) throw std::runtime_error("cannot write output CSV");
  out << std::setprecision(17);
  out << "# source=redesigned Albert_Young case, 80x80, 10 steps, dt=0.002 days\n";
  out << "# n1=" << n1 << "\n# n2=" << n2 << "\n# steps=" << steps
      << "\n# dt=" << dt << "\n";
  out << "# min0=" << *std::min_element(f0.begin(), f0.end()) << "\n";
  out << "# min10=" << *std::min_element(f10.begin(), f10.end()) << "\n";
  out << "i,j,f0,f10\n";
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      const std::size_t k = i * n2 + j;
      out << i << "," << j << "," << f0[k] << "," << f10[k] << "\n";
    }
  }
  std::cout << "min0=" << *std::min_element(f0.begin(), f0.end())
            << " min10=" << *std::min_element(f10.begin(), f10.end()) << "\n";
  return 0;
}
