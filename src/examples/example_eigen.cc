// Minimal embedding example using Eigen::MatrixXd directly.
//
// This is the complete, compilable version of the embedding snippet in the
// project README. It uses the direct Eigen coefficient and field interfaces.

#include <cstddef>
#include <iostream>

#include <Eigen/Dense>

#include "boundary_conditions.h"
#include "equation.h"
#include "mesh.h"
#include "solver.h"

int main() {
  using namespace sayram2d;

  Mesh mesh({0.0, 0.3, 0.7, 1.0}, {0.0, 0.4, 1.0, 1.7});

  const std::size_t n1 = mesh.n1();
  const std::size_t n2 = mesh.n2();

  Eigen::MatrixXd G = Eigen::MatrixXd::Constant(n1, n2, 1.0);
  Eigen::MatrixXd D11 = Eigen::MatrixXd::Constant(n1, n2, 1.0);
  Eigen::MatrixXd D12 = Eigen::MatrixXd::Zero(n1, n2);
  Eigen::MatrixXd D22 = Eigen::MatrixXd::Constant(n1, n2, 1.0);

  BoundaryConditions bc;
  bc.set_dirichlet(Side::Q1Min, 0.5);
  bc.set_zero_flux(Side::Q1Max);
  bc.set_zero_flux(Side::Q2Min);
  bc.set_zero_flux(Side::Q2Max);

  Equation equation(G, D11, D12, D22, bc);
  // Optional; inverse loss defaults to zero.
  equation.set_loss(0.0);

  Solver solver(mesh);
  const Eigen::MatrixXd f_in = Eigen::MatrixXd::Ones(n1, n2);
  Eigen::MatrixXd f_out(n1, n2);
  const double dt = 0.01;
  solver.solve(equation, f_in, f_out, dt);

  std::cout << "example_eigen: dt=" << dt << " min f=" << f_out.minCoeff()
            << " max f=" << f_out.maxCoeff() << "\n";
  return 0;
}
