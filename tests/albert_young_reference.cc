// Run the real Albert-Young case on the preprocessed diffusion data.
// Usage: albert_young_reference <converted.bin> <output.csv>

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

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

  const Eigen::MatrixXd f0 = data.f();
  Eigen::MatrixXd f_in = f0;
  Eigen::MatrixXd f_out(n1, n2);
  for (int step = 0; step < steps; ++step) {
    solver.solve(data.equation(), f_in, f_out, dt);
    std::swap(f_in, f_out);
  }
  const Eigen::MatrixXd& f10 = f_in;

  std::ofstream out(argv[2]);
  if (!out) throw std::runtime_error("cannot write output CSV");
  out << std::setprecision(17);
  out << "# source=redesigned Albert_Young case, " << n1 << "x" << n2
      << ", " << steps << " steps, dt=" << dt << " days\n";
  out << "# n1=" << n1 << "\n# n2=" << n2 << "\n# steps=" << steps
      << "\n# dt=" << dt << "\n";
  out << "# index_order=row-major i*n2+j; i=q1(alpha0), j=q2(logE)\n";
  out << "# q1_unit=rad, q2=log(E/MeV)\n";
  out << "# min0=" << f0.minCoeff() << "\n";
  out << "# min10=" << f10.minCoeff() << "\n";
  out << "i,j,f0,f10\n";
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      out << i << "," << j << "," << f0(static_cast<Eigen::Index>(i),
                                        static_cast<Eigen::Index>(j))
          << "," << f10(static_cast<Eigen::Index>(i),
                         static_cast<Eigen::Index>(j))
          << "\n";
    }
  }
  std::cout << "min0=" << f0.minCoeff() << " min10=" << f10.minCoeff()
            << "\n";
  return 0;
}
