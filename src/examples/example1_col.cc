// Example 1 with external column-major storage.
//
// This is the same PPFV test problem as example1.cc. The coefficients and
// fields are flat std::vector<double> buffers instead of Array2D. The host
// index is column-major:
//
//   data[i + j*n1]
//
// so the ArrayView2D strides are (1, n1). The views are constructed first;
// Equation and Solver then see only the logical (q1, q2) field.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

#include "array_view2d.h"
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
  const std::size_t n_cells = n1 * n2;

  std::vector<double> G(n_cells, 1.0);
  std::vector<double> D11(n_cells, 1.0);
  std::vector<double> D12(n_cells, 0.8);
  std::vector<double> D22(n_cells, 1.0);

  // Column-major host layout and matching ArrayView2D strides.
  const std::ptrdiff_t stride1 = 1;
  const std::ptrdiff_t stride2 = static_cast<std::ptrdiff_t>(n1);

  // Construct every view before Equation/Solver.
  auto G_view = make_view2d(G.data(), n1, n2, stride1, stride2);
  auto D11_view = make_view2d(D11.data(), n1, n2, stride1, stride2);
  auto D12_view = make_view2d(D12.data(), n1, n2, stride1, stride2);
  auto D22_view = make_view2d(D22.data(), n1, n2, stride1, stride2);

  std::vector<double> f_in(n_cells);
  std::vector<double> f_out(n_cells, 0.0);
  auto f_in_view = make_view2d(f_in.data(), n1, n2, stride1, stride2);
  auto f_out_view = make_view2d(f_out.data(), n1, n2, stride1, stride2);

  for (std::size_t j = 0; j < n2; ++j) {
    for (std::size_t i = 0; i < n1; ++i) {
      const double x = mesh.q1_center(i);
      const double y = mesh.q2_center(j);
      f_in[i + j * n1] =
          std::exp(-y / 0.02) * std::sin(0.5 * physics::pi * x);
    }
  }

  std::vector<double> q2_min_values(n1 + 1);
  for (std::size_t i = 0; i <= n1; ++i) {
    q2_min_values[i] = std::sin(0.5 * physics::pi * mesh.q1_edge(i));
  }

  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.0);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_dirichlet(Side::Q2Min, q2_min_values);
  bc.set_dirichlet(Side::Q2Max, 0.0);

  Equation equation(G_view, D11_view, D12_view, D22_view, bc);
  Solver solver(mesh);

  const double dt = 1.0e-4;
  const int steps = 100;  // final time = 0.01
  for (int step = 1; step <= steps; ++step) {
    solver.solve(equation, f_in_view, f_out_view, dt);
    // Swap the views, not the vectors, so both views stay valid.
    std::swap(f_in_view, f_out_view);
    if (step == 20 || step == steps) {
      double min = f_in_view(0, 0);
      double max = min;
      for (std::size_t j = 0; j < n2; ++j) {
        for (std::size_t i = 0; i < n1; ++i) {
          min = std::min(min, f_in_view(i, j));
          max = std::max(max, f_in_view(i, j));
        }
      }
      std::cout << "column-major: t = " << step * dt
                << "  min f = " << min << "  max f = " << max << "\n";
    }
  }

  return 0;
}
