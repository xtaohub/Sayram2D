// Time-dependent coefficients and boundary conditions.
//
// This example shows that an Equation can be constructed once while its
// coefficient views and BoundaryConditions reference stay live. The diffusion
// coefficients and the Q1Min Dirichlet array are updated before each solve.

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
  const double dt = 0.01;
  const int steps = 100;       // final time = 1.0
  const double tau = 10.0;     // slow decay: exp(-1/10) ~= 0.905 at t=1

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
  Array2D<double> D11(n1, n2);
  Array2D<double> D12(n1, n2);
  Array2D<double> D22(n1, n2);
  Array2D<double> f_in(n1, n2, 1.0);
  Array2D<double> f_out(n1, n2, 0.0);

  // Q1Min values live at the n2+1 vertices along the boundary.
  std::vector<double> q1_min_values(n2 + 1);
  for (std::size_t j = 0; j <= n2; ++j) {
    q1_min_values[j] =
        0.5 * (1.0 + 0.25 * std::sin(physics::pi * mesh.q2_edge(j)));
  }

  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, q1_min_values);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_zero_flux(Side::Q2Min);
  bc.set_zero_flux(Side::Q2Max);

  Equation equation(G.view(), D11.view(), D12.view(), D22.view(), bc);
  Solver solver(mesh);

  auto update_parameters = [&](double time) {
    const double decay = std::exp(-time / tau);
    D11.fill(decay);
    D12.fill(0.8 * decay);
    D22.fill(decay);
    for (std::size_t j = 0; j <= n2; ++j) {
      q1_min_values[j] =
          0.5 * decay *
          (1.0 + 0.25 * std::sin(physics::pi * mesh.q2_edge(j)));
    }
    bc.set_dirichlet(Side::Q1Min, q1_min_values);
  };

  auto print_state = [&](double time) {
    const auto minmax = std::minmax_element(f_in.begin(), f_in.end());
    std::cout << "t=" << time << " D11=" << D11(0, 0)
              << " g_D(mid)=" << bc.value(Side::Q1Min, n2 / 2)
              << " min f=" << *minmax.first
              << " max f=" << *minmax.second << "\n";
  };

  update_parameters(0.0);
  print_state(0.0);

  for (int step = 1; step <= steps; ++step) {
    const double time = static_cast<double>(step) * dt;
    // The Equation sees these updates because it holds non-owning views and a
    // reference to bc.
    update_parameters(time);
    solver.solve(equation, f_in.view(), f_out.view(), dt);
    std::swap(f_in, f_out);
    if (step == steps / 2 || step == steps) {
      print_state(time);
    }
  }

  return 0;
}
