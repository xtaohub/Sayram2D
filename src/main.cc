// A short introduction to the Sayram2D library.
//
// It solves a small 2-D diffusion problem with constant diagonal
// coefficients and mixed Dirichlet/ZeroFlux boundaries. The same pattern is
// used for realistic cases: build a Mesh, describe the coefficient arrays and
// boundary conditions, construct Equation and Solver, then call solve() once
// per time step.

#include <algorithm>
#include <iostream>
#include <vector>

#include "array_view2d.h"
#include "boundary_conditions.h"
#include "equation.h"
#include "mesh.h"
#include "solver.h"

int main() {
  using namespace sayram2d;

  // 4 x 2 orthogonal grid; spacing may be nonuniform.
  Mesh mesh({0.0, 0.25, 0.5, 0.75, 1.0}, {0.0, 0.5, 1.0});

  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();
  const std::size_t cells = n1 * n2;

  // Cell-centered coefficients. D12 is the symmetric off-diagonal component.
  std::vector<double> G(cells, 1.0);
  std::vector<double> D11(cells, 1.0);
  std::vector<double> D12(cells, 0.0);
  std::vector<double> D22(cells, 1.0);

  // f = 0 on the lower-q1 side; all other sides are zero flux.
  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.0);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_zero_flux(Side::Q2Min);
  bc.set_zero_flux(Side::Q2Max);

  Equation equation(make_view2d(G, n1, n2), make_view2d(D11, n1, n2),
                    make_view2d(D12, n1, n2), make_view2d(D22, n1, n2), bc);

  // The caller owns the evolving field. f_in is time t; f_out receives t+dt.
  std::vector<double> f_in(cells, 1.0);
  std::vector<double> f_out(cells, 0.0);

  Solver solver(mesh);
  const double dt = 0.01;
  for (int step = 0; step < 5; ++step) {
    solver.solve(equation, make_view2d(f_in, n1, n2),
                 make_view2d(f_out, n1, n2), dt);
    std::swap(f_in, f_out);
  }

  const auto minmax = std::minmax_element(f_in.begin(), f_in.end());
  const std::size_t center = (n1 / 2) * n2 + (n2 / 2);
  std::cout << "center f = " << f_in[center] << "\n"
            << "min f    = " << *minmax.first << "\n"
            << "max f    = " << *minmax.second << "\n";
  return 0;
}
