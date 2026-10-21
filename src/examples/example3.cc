// Time-dependent coefficients and boundary conditions.
//
// This example shows that an Equation can be constructed once while its
// coefficient views and BoundaryConditions reference stay live. The diffusion
// coefficients and the Q1Min Dirichlet value are updated in place before each
// solve.

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
  const std::size_t cells = n1 * n2;
  const Shape2 shape{n1, n2};

  std::vector<double> G(cells, 1.0);
  std::vector<double> D11(cells);
  std::vector<double> D12(cells);
  std::vector<double> D22(cells);
  std::vector<double> f_in(cells, 1.0);
  std::vector<double> f_out(cells, 0.0);

  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.5);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_zero_flux(Side::Q2Min);
  bc.set_zero_flux(Side::Q2Max);

  Equation equation(make_view2d(G, shape), make_view2d(D11, shape),
                    make_view2d(D12, shape), make_view2d(D22, shape), bc);
  Solver solver(mesh);

  auto update_parameters = [&](double time) {
    const double decay = std::exp(-time / tau);
    std::fill(D11.begin(), D11.end(), decay);
    std::fill(D12.begin(), D12.end(), 0.25 * decay);
    std::fill(D22.begin(), D22.end(), decay);
    bc.set_dirichlet(Side::Q1Min, 0.5 * decay);
  };

  auto print_state = [&](double time) {
    const auto minmax = std::minmax_element(f_in.begin(), f_in.end());
    std::cout << "t=" << time << " D11=" << D11[0]
              << " g_D=" << bc.value(Side::Q1Min, 0)
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
    solver.solve(equation, make_view2d(f_in, shape),
                 make_view2d(f_out, shape), dt);
    std::swap(f_in, f_out);
    if (step == steps / 2 || step == steps) {
      print_state(time);
    }
  }

  return 0;
}
