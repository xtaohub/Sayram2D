#include "sayram2d/solver.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace sayram2d {
namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();

double dot(Vec2 a, Vec2 b) noexcept { return a.q1 * b.q1 + a.q2 * b.q2; }

}  // namespace

Solver::Solver(const Mesh& mesh)
    : mesh_(mesh),
      matrix_(static_cast<int>(mesh.n_cells()), static_cast<int>(mesh.n_cells())),
      rhs_(static_cast<int>(mesh.n_cells())),
      solution_(static_cast<int>(mesh.n_cells())),
      vertex_f_((mesh.n1() + 1) * (mesh.n2() + 1)),
      lambda_(mesh.n_cells()) {
  rhs_.setZero();
  solution_.setZero();
  triplets_.reserve(6 * mesh.n_cells());
  build_pattern();
}

void Solver::build_pattern() {
  triplets_.clear();
  triplets_.reserve(6 * mesh_.n_cells());

  // Every cell has a diagonal entry.
  for (std::size_t i = 0; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      const auto k = static_cast<int>(mesh_.cell_index({i, j}));
      triplets_.emplace_back(k, k, 1.0);
    }
  }

  // Interior-face pattern: four entries per shared face. Boundary Dirichlet
  // faces add only a diagonal entry, which is already present.
  for (std::size_t i = 1; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      const auto k = static_cast<int>(mesh_.cell_index({i, j}));
      const auto l = static_cast<int>(mesh_.cell_index({i - 1, j}));
      triplets_.emplace_back(k, k, 1.0);
      triplets_.emplace_back(k, l, 1.0);
      triplets_.emplace_back(l, l, 1.0);
      triplets_.emplace_back(l, k, 1.0);
    }
  }
  for (std::size_t j = 1; j < mesh_.n2(); ++j) {
    for (std::size_t i = 0; i < mesh_.n1(); ++i) {
      const auto k = static_cast<int>(mesh_.cell_index({i, j}));
      const auto l = static_cast<int>(mesh_.cell_index({i, j - 1}));
      triplets_.emplace_back(k, k, 1.0);
      triplets_.emplace_back(k, l, 1.0);
      triplets_.emplace_back(l, l, 1.0);
      triplets_.emplace_back(l, k, 1.0);
    }
  }

  matrix_.setFromTriplets(triplets_.begin(), triplets_.end());
  matrix_.makeCompressed();
  linear_solver_.analyzePattern(matrix_);
  triplets_.clear();
}

void Solver::validate(const Equation& equation, ArrayView2D<const double> f_in,
                      ArrayView2D<double> f_out, double dt) const {
  if (!(dt > 0.0) || !std::isfinite(dt)) {
    throw std::invalid_argument("Solver::solve: dt must be positive and finite");
  }
  if (f_in.extent1() != mesh_.n1() || f_in.extent2() != mesh_.n2()) {
    throw std::invalid_argument("Solver::solve: f_in has the wrong extents");
  }
  if (f_out.extent1() != mesh_.n1() || f_out.extent2() != mesh_.n2()) {
    throw std::invalid_argument("Solver::solve: f_out has the wrong extents");
  }
  if (f_in.data() == nullptr || f_out.data() == nullptr) {
    throw std::invalid_argument("Solver::solve: input/output view is null");
  }

  const auto G = equation.G_view();
  const auto D11 = equation.D11_view();
  const auto D12 = equation.D12_view();
  const auto D22 = equation.D22_view();
  if (G.extent1() != mesh_.n1() || G.extent2() != mesh_.n2() ||
      D11.extent1() != mesh_.n1() || D11.extent2() != mesh_.n2() ||
      D12.extent1() != mesh_.n1() || D12.extent2() != mesh_.n2() ||
      D22.extent1() != mesh_.n1() || D22.extent2() != mesh_.n2()) {
    throw std::invalid_argument(
        "Solver::solve: coefficient views do not match the mesh");
  }
  if (G.data() == nullptr || D11.data() == nullptr || D12.data() == nullptr ||
      D22.data() == nullptr) {
    throw std::invalid_argument("Solver::solve: coefficient view is null");
  }

  if (equation.has_loss() && equation.loss_is_constant() &&
      !std::isfinite(equation.constant_inv_tau())) {
    throw std::invalid_argument("Solver::solve: inverse loss is not finite");
  }
  if (equation.has_loss() && !equation.loss_is_constant()) {
    const auto tau = equation.inv_tau_view();
    if (tau.extent1() != mesh_.n1() || tau.extent2() != mesh_.n2() ||
        tau.data() == nullptr) {
      throw std::invalid_argument(
          "Solver::solve: inverse loss view does not match the mesh");
    }
  }

  for (std::size_t i = 0; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      const double g = G(i, j);
      const double d11 = D11(i, j);
      const double d12 = D12(i, j);
      const double d22 = D22(i, j);
      const double f = f_in(i, j);
      const double tau = equation.inv_tau(i, j);
      if (!std::isfinite(g) || !std::isfinite(d11) || !std::isfinite(d12) ||
          !std::isfinite(d22) || !std::isfinite(f) || !std::isfinite(tau)) {
        throw std::invalid_argument(
            "Solver::solve: non-finite coefficient, input, or loss value");
      }
      if (!(g > 0.0)) {
        throw std::invalid_argument("Solver::solve: G must be positive");
      }
      if (!(d11 > 0.0) || !(d22 > 0.0) ||
          !(d11 * d22 - d12 * d12 > 0.0)) {
        throw std::invalid_argument(
            "Solver::solve: diffusion tensor must be positive definite");
      }
    }
  }

  equation.bc().validate(mesh_);
}

void Solver::update_lambda(const Equation& equation) {
  for (std::size_t i = 0; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      const double g = equation.G(i, j);
      const double d11 = equation.D11(i, j);
      const double d12 = equation.D12(i, j);
      const double d22 = equation.D22(i, j);
      lambda_[mesh_.cell_index({i, j})] = {g * d11, g * d12, g * d12,
                                           g * d22};
    }
  }
}

void Solver::update_vertex_f(const Equation& equation,
                             ArrayView2D<const double> f_in) {
  const std::size_t n1 = mesh_.n1();
  const std::size_t n2 = mesh_.n2();
  const std::size_t nv1 = n1 + 1;
  auto v = [&](std::size_t i, std::size_t j) -> double& {
    return vertex_f_[j * nv1 + i];
  };
  auto lerp = [](double x, double x0, double x1, double f0, double f1) {
    const double denom = x1 - x0;
    if (denom == 0.0) return 0.5 * (f0 + f1);
    const double w0 = (x1 - x) / denom;
    const double w1 = (x - x0) / denom;
    return w0 * f0 + w1 * f1;
  };

  // Corners: the single adjacent cell value.
  v(0, 0) = f_in(0, 0);
  v(n1, 0) = f_in(n1 - 1, 0);
  v(0, n2) = f_in(0, n2 - 1);
  v(n1, n2) = f_in(n1 - 1, n2 - 1);

  // Interior vertices: bilinear interpolation from four surrounding cells.
  for (std::size_t i = 1; i < n1; ++i) {
    for (std::size_t j = 1; j < n2; ++j) {
      const double q1v = mesh_.q1_edge(i);
      const double q2v = mesh_.q2_edge(j);
      const double q1l = mesh_.q1_center(i - 1);
      const double q1r = mesh_.q1_center(i);
      const double q2l = mesh_.q2_center(j - 1);
      const double q2r = mesh_.q2_center(j);
      const double wx_l = (q1r - q1v) / (q1r - q1l);
      const double wx_r = (q1v - q1l) / (q1r - q1l);
      const double wy_l = (q2r - q2v) / (q2r - q2l);
      const double wy_r = (q2v - q2l) / (q2r - q2l);
      v(i, j) = wx_l * wy_l * f_in(i - 1, j - 1) +
                wx_r * wy_l * f_in(i, j - 1) +
                wx_l * wy_r * f_in(i - 1, j) +
                wx_r * wy_r * f_in(i, j);
    }
  }

  // Boundary edges parallel to q2 (Q1Min/Q1Max), excluding corners.
  for (std::size_t j = 1; j < n2; ++j) {
    const double q2v = mesh_.q2_edge(j);
    const double q2l = mesh_.q2_center(j - 1);
    const double q2r = mesh_.q2_center(j);
    v(0, j) = lerp(q2v, q2l, q2r, f_in(0, j - 1), f_in(0, j));
    v(n1, j) = lerp(q2v, q2l, q2r, f_in(n1 - 1, j - 1),
                    f_in(n1 - 1, j));
  }

  // Boundary edges parallel to q1 (Q2Min/Q2Max), excluding corners.
  for (std::size_t i = 1; i < n1; ++i) {
    const double q1v = mesh_.q1_edge(i);
    const double q1l = mesh_.q1_center(i - 1);
    const double q1r = mesh_.q1_center(i);
    v(i, 0) = lerp(q1v, q1l, q1r, f_in(i - 1, 0), f_in(i, 0));
    v(i, n2) = lerp(q1v, q1l, q1r, f_in(i - 1, n2 - 1),
                    f_in(i, n2 - 1));
  }

  // Dirichlet data overwrite the boundary vertices. Corner consistency was
  // checked by BoundaryConditions::validate before this function is called.
  const BoundaryConditions& bc = equation.bc();
  for (Side side : all_sides()) {
    if (!bc.is_dirichlet(side)) continue;
    switch (side) {
      case Side::Q1Min:
        for (std::size_t j = 0; j <= n2; ++j) v(0, j) = bc.value(side, j);
        break;
      case Side::Q1Max:
        for (std::size_t j = 0; j <= n2; ++j) v(n1, j) = bc.value(side, j);
        break;
      case Side::Q2Min:
        for (std::size_t i = 0; i <= n1; ++i) v(i, 0) = bc.value(side, i);
        break;
      case Side::Q2Max:
        for (std::size_t i = 0; i <= n1; ++i) v(i, n2) = bc.value(side, i);
        break;
    }
  }
}

double Solver::alpha_one_sided(const Face& face, const Mat2& lambda,
                               double* alpha_a, double* alpha_b) const {
  const Vec2 center = mesh_.cell_center(face.owner);
  const Vec2 vka{face.a.q1 - center.q1, face.a.q2 - center.q2};
  const Vec2 vkb{face.b.q1 - center.q1, face.b.q2 - center.q2};
  const Vec2 rvka{vka.q2, -vka.q1};
  const Vec2 rvkb{vkb.q2, -vkb.q1};

  const auto lambda_times = [&](Vec2 r) {
    return Vec2{lambda.a00 * r.q1 + lambda.a01 * r.q2,
                lambda.a10 * r.q1 + lambda.a11 * r.q2};
  };

  const double numerator_a = dot(face.normal, lambda_times(rvkb));
  const double numerator_b = dot(face.normal, lambda_times(rvka));
  const double denominator_a = dot(vka, rvkb);
  const double denominator_b = dot(vkb, rvka);

  *alpha_a = face.length * numerator_a / denominator_a;
  *alpha_b = face.length * numerator_b / denominator_b;
  return (*alpha_a) * vertex_value(face.va) +
         (*alpha_b) * vertex_value(face.vb);
}

void Solver::add_face_pair(const Face& face_k, const Face& face_l,
                           ArrayView2D<const double> f_in) {
  const Ind k = face_k.owner;
  const Ind l = face_l.owner;
  const Mat2& lambda_k = lambda_[mesh_.cell_index(k)];
  const Mat2& lambda_l = lambda_[mesh_.cell_index(l)];

  double alpha_ka = 0.0;
  double alpha_kb = 0.0;
  double alpha_la = 0.0;
  double alpha_lb = 0.0;
  const double a_sigma_k =
      alpha_one_sided(face_k, lambda_k, &alpha_ka, &alpha_kb);
  const double a_sigma_l =
      alpha_one_sided(face_l, lambda_l, &alpha_la, &alpha_lb);

  const double denominator =
      std::abs(a_sigma_k) + std::abs(a_sigma_l) + 2.0 * kEps;
  const double mu_k = (std::abs(a_sigma_l) + kEps) / denominator;
  const double mu_l = 1.0 - mu_k;

  const double b_sigma = mu_l * a_sigma_l - mu_k * a_sigma_k;
  const double b_abs = std::abs(b_sigma);
  const double b_plus = (b_abs + b_sigma) / 2.0;
  const double b_minus = (b_abs - b_sigma) / 2.0;

  const double f_k = f_in(k.i, k.j);
  const double f_l = f_in(l.i, l.j);
  const double A_k =
      mu_k * (alpha_ka + alpha_kb) + b_plus / (f_k + kEps);
  const double A_l =
      mu_l * (alpha_la + alpha_lb) + b_minus / (f_l + kEps);

  const auto K = static_cast<int>(mesh_.cell_index(k));
  const auto L = static_cast<int>(mesh_.cell_index(l));
  triplets_.emplace_back(K, K, A_k);
  triplets_.emplace_back(K, L, -A_l);
  triplets_.emplace_back(L, L, A_l);
  triplets_.emplace_back(L, K, -A_k);
}

void Solver::assemble(const Equation& equation, ArrayView2D<const double> f_in,
                      double dt) {
  rhs_.setZero();
  triplets_.clear();

  // Each interior face is processed once. Q1Min faces for i >= 1 and Q2Min
  // faces for j >= 1 cover every interior face exactly once.
  for (std::size_t i = 1; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      const Ind k{i, j};
      const Face face_k = mesh_.face(k, Side::Q1Min);
      assert(!face_k.boundary);
      const Face face_l = mesh_.face(face_k.neighbor, Side::Q1Max);
      add_face_pair(face_k, face_l, f_in);
    }
  }
  for (std::size_t j = 1; j < mesh_.n2(); ++j) {
    for (std::size_t i = 0; i < mesh_.n1(); ++i) {
      const Ind k{i, j};
      const Face face_k = mesh_.face(k, Side::Q2Min);
      assert(!face_k.boundary);
      const Face face_l = mesh_.face(face_k.neighbor, Side::Q2Max);
      add_face_pair(face_k, face_l, f_in);
    }
  }

  // Dirichlet boundary faces. ZeroFlux contributes no face term.
  const BoundaryConditions& bc = equation.bc();
  for (Side side : all_sides()) {
    if (!bc.is_dirichlet(side)) continue;

    std::size_t count = 0;
    switch (side) {
      case Side::Q1Min:
      case Side::Q1Max:
        count = mesh_.n2();
        break;
      case Side::Q2Min:
      case Side::Q2Max:
        count = mesh_.n1();
        break;
    }

    for (std::size_t k = 0; k < count; ++k) {
      Ind cell{};
      switch (side) {
        case Side::Q1Min:
          cell = {0, k};
          break;
        case Side::Q1Max:
          cell = {mesh_.n1() - 1, k};
          break;
        case Side::Q2Min:
          cell = {k, 0};
          break;
        case Side::Q2Max:
          cell = {k, mesh_.n2() - 1};
          break;
      }

      const Face face = mesh_.face(cell, side);
      const Mat2& lambda = lambda_[mesh_.cell_index(cell)];
      double alpha_a = 0.0;
      double alpha_b = 0.0;
      const double a_sigma =
          alpha_one_sided(face, lambda, &alpha_a, &alpha_b);
      const double b_sigma = -a_sigma;
      const double b_abs = std::abs(b_sigma);
      const double b_plus = (b_abs + b_sigma) / 2.0;
      const double b_minus = (b_abs - b_sigma) / 2.0;
      const double f_k = f_in(cell.i, cell.j);
      const double A_k = alpha_a + alpha_b + b_plus / (f_k + kEps);
      const auto K = static_cast<int>(mesh_.cell_index(cell));
      triplets_.emplace_back(K, K, A_k);
      rhs_(K) += b_minus;
    }
  }

  // Time term, implicit loss (tau > 0), and explicit source (tau < 0).
  for (std::size_t i = 0; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      const Ind cell{i, j};
      const double uk = equation.G(i, j) * mesh_.cell_measure(cell) / dt;
      const double f_old = f_in(i, j);
      const double tau = equation.inv_tau(i, j);
      const auto K = static_cast<int>(mesh_.cell_index(cell));
      if (tau >= 0.0) {
        triplets_.emplace_back(K, K, uk * (1.0 + dt * tau));
        rhs_(K) += uk * f_old;
      } else {
        // -f/tau is a source: R += G*area*f/|tau| = uk*dt*(-tau)*f.
        triplets_.emplace_back(K, K, uk);
        rhs_(K) += uk * f_old * (1.0 - dt * tau);
      }
    }
  }

  matrix_.setFromTriplets(triplets_.begin(), triplets_.end());
  matrix_.makeCompressed();
}

void Solver::solve(const Equation& equation, ArrayView2D<const double> f_in,
                   ArrayView2D<double> f_out, double dt) {
  validate(equation, f_in, f_out, dt);

  update_lambda(equation);
  update_vertex_f(equation, f_in);
  assemble(equation, f_in, dt);

  linear_solver_.factorize(matrix_);
  if (linear_solver_.info() != Eigen::Success) {
    throw std::runtime_error("Solver::solve: sparse factorization failed");
  }
  solution_ = linear_solver_.solve(rhs_);
  if (linear_solver_.info() != Eigen::Success) {
    throw std::runtime_error("Solver::solve: sparse solve failed");
  }

  // All reads of f_in are complete at this point, so f_out may alias f_in.
  for (std::size_t i = 0; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      f_out(i, j) = solution_(static_cast<int>(mesh_.cell_index({i, j})));
    }
  }
}

}  // namespace sayram2d
