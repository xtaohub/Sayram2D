#pragma once

#include <cstddef>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "array_view2d.h"
#include "equation.h"
#include "mesh.h"

namespace sayram2d {

class Solver {
 public:
  explicit Solver(const Mesh& mesh);

  void solve(const Equation& equation, ArrayView2D<const double> f_in,
             ArrayView2D<double> f_out, double dt);

 private:
  struct Mat2 {
    double a00 = 0.0;
    double a01 = 0.0;
    double a10 = 0.0;
    double a11 = 0.0;
  };

  using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor>;
  using Triplet = Eigen::Triplet<double>;
  using LinearSolver =
      Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>;

  const Mesh& mesh_;
  SparseMatrix matrix_;
  LinearSolver linear_solver_;
  std::vector<Triplet> triplets_;
  Eigen::VectorXd rhs_;
  Eigen::VectorXd solution_;

  // Vertex values at (n1+1) x (n2+1) grid vertices.
  std::vector<double> vertex_f_;
  // Cell-centered Lambda = G * D, one symmetric 2x2 tensor per cell.
  std::vector<Mat2> lambda_;

  void build_pattern();

  void validate(const Equation& equation, ArrayView2D<const double> f_in,
                ArrayView2D<double> f_out, double dt) const;
  void update_lambda(const Equation& equation);
  void update_vertex_f(const Equation& equation,
                       ArrayView2D<const double> f_in);
  void assemble(const Equation& equation, ArrayView2D<const double> f_in,
                double dt);
  void add_face_pair(const Face& face_k, const Face& face_l,
                     ArrayView2D<const double> f_in);
  double alpha_one_sided(const Face& face, const Mat2& lambda,
                         double* alpha_a, double* alpha_b) const;
  double vertex_value(VtxInd vertex) const noexcept {
    return vertex_f_[flat_index(vertex.i, vertex.j, mesh_.n2() + 1)];
  }
};

}  // namespace sayram2d
