// A short introduction to the Sayram2D library.
//
// It solves a small 2-D diffusion problem with constant diagonal
// coefficients and mixed Dirichlet/ZeroFlux boundaries. The same pattern is
// used for realistic cases: build a Mesh, provide Eigen coefficient matrices
// and boundary conditions, construct Equation and Solver, then call solve()
// once per time step.

#include <Eigen/Dense>

#include <iostream>
#include <utility>

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

  // Cell-centered coefficients. Eigen::MatrixXd is column-major and is
  // accepted directly by Equation without a copy.
  Eigen::MatrixXd G = Eigen::MatrixXd::Constant(n1, n2, 1.0);
  Eigen::MatrixXd D11 = Eigen::MatrixXd::Constant(n1, n2, 1.0);
  Eigen::MatrixXd D12 = Eigen::MatrixXd::Zero(n1, n2);
  Eigen::MatrixXd D22 = Eigen::MatrixXd::Constant(n1, n2, 1.0);

  // f = 0 on the lower-q1 side; all other sides are zero flux.
  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.0);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_zero_flux(Side::Q2Min);
  bc.set_zero_flux(Side::Q2Max);

  Equation equation(G, D11, D12, D22, bc);

  // The caller owns the evolving field. f_in is time t; f_out receives t+dt.
  Eigen::MatrixXd f_in = Eigen::MatrixXd::Ones(n1, n2);
  Eigen::MatrixXd f_out = Eigen::MatrixXd::Zero(n1, n2);

  Solver solver(mesh);
  const double dt = 0.01;
  for (int step = 0; step < 5; ++step) {
    solver.solve(equation, f_in, f_out, dt);
    std::swap(f_in, f_out);
  }

  std::cout << "center f = " << f_in(n1 / 2, n2 / 2) << "\n"
            << "min f    = " << f_in.minCoeff() << "\n"
            << "max f    = " << f_in.maxCoeff() << "\n";
  return 0;
}
