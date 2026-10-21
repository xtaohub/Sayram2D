// Time-independent test problem from the PPFV paper.
//
//   f_t = div(D grad f),  (x,y) in [0,1]^2,
//   f(x,y,0) = exp(-y/0.02) sin(pi x / 2),
//   f(0,y,t) = 0,  f_x(1,y,t) = 0,
//   f(x,0,t) = f(x,0,0),  f(x,1,t) = 0,
//
// with G = D11 = D22 = 1 and D12 = 0.8. The example uses Array2D for owned
// coefficient and field storage and ArrayView2D through Equation/Solver.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

#include "array2d.h"
#include "boundary_conditions.h"
#include "equation.h"
#include "mesh.h"
#include "physics.h"
#include "solver.h"

int main() {
  using namespace sayram2d;

  const std::size_t n = 50;

  std::vector<double> q1_edges(n + 1);
  std::vector<double> q2_edges(n + 1);
  for (std::size_t i = 0; i <= n; ++i) {
    q1_edges[i] = static_cast<double>(i) / static_cast<double>(n);
    q2_edges[i] = static_cast<double>(i) / static_cast<double>(n);
  }
  Mesh mesh(std::move(q1_edges), std::move(q2_edges));

  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();

  Array2D<double> G(n1, n2, 1.0);
  Array2D<double> D11(n1, n2, 1.0);
  Array2D<double> D12(n1, n2, 0.8);
  Array2D<double> D22(n1, n2, 1.0);

  std::vector<double> q2_min_values(n1 + 1);
  for (std::size_t i = 0; i <= n1; ++i) {
    q2_min_values[i] = std::sin(0.5 * physics::pi * mesh.q1_edge(i));
  }

  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.0);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_dirichlet(Side::Q2Min, q2_min_values);
  bc.set_dirichlet(Side::Q2Max, 0.0);

  Equation equation(G.view(), D11.view(), D12.view(), D22.view(), bc);
  Solver solver(mesh);

  Array2D<double> f_in(n1, n2);
  Array2D<double> f_out(n1, n2, 0.0);
  for (std::size_t j = 0; j < n2; ++j) {
    for (std::size_t i = 0; i < n1; ++i) {
      const double x = mesh.q1_center(i);
      const double y = mesh.q2_center(j);
      f_in(i, j) = std::exp(-y / 0.02) * std::sin(0.5 * physics::pi * x);
    }
  }

  const double dt = 1.0e-4;
  const int steps = 100;  // final time = 0.01
  for (int step = 1; step <= steps; ++step) {
    solver.solve(equation, f_in.view(), f_out.view(), dt);
    std::swap(f_in, f_out);
    if (step == 20 || step == steps) {
      const auto minmax = std::minmax_element(f_in.begin(), f_in.end());
      std::cout << "t = " << step * dt << "  min f = " << *minmax.first
                << "  max f = " << *minmax.second << "\n";
    }
  }

  return 0;
}
