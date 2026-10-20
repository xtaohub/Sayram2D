#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

#include "array_view2d.h"
#include "boundary_conditions.h"
#include "equation.h"
#include "mesh.h"
#include "solver.h"

namespace reference {

constexpr double pi = 3.141592653589793238462643383279502884;

struct CaseConfig {
  std::array<sayram2d::BCType, 4> bc{};
  double inv_tau = 0.0;
  bool off_diagonal = false;
  bool varying_coefficients = true;
};

inline std::array<sayram2d::BCType, 4> all_dirichlet() {
  using sayram2d::BCType;
  return {BCType::Dirichlet, BCType::Dirichlet, BCType::Dirichlet,
          BCType::Dirichlet};
}

inline std::array<sayram2d::BCType, 4> mixed_dirichlet_zeroflux() {
  using sayram2d::BCType;
  return {BCType::Dirichlet, BCType::ZeroFlux, BCType::Dirichlet,
          BCType::ZeroFlux};
}

inline std::array<sayram2d::BCType, 4> all_zeroflux() {
  using sayram2d::BCType;
  return {BCType::ZeroFlux, BCType::ZeroFlux, BCType::ZeroFlux,
          BCType::ZeroFlux};
}

inline std::vector<double> uniform_edges(double a, double b, std::size_t n) {
  if (n == 0) throw std::runtime_error("uniform_edges: n must be positive");
  std::vector<double> edges(n + 1);
  const double h = (b - a) / static_cast<double>(n);
  for (std::size_t k = 0; k <= n; ++k) {
    edges[k] = a + h * static_cast<double>(k);
  }
  return edges;
}

inline std::vector<double> nonuniform_q1_edges() {
  return {0.00, 0.03, 0.08, 0.18, 0.32, 0.50, 0.70, 0.85, 0.95, 0.99, 1.00};
}

inline std::vector<double> nonuniform_q2_edges() {
  return {0.00, 0.05, 0.12, 0.22, 0.35, 0.50, 0.65, 0.78, 0.88, 0.95, 1.00,
          1.20};
}

inline double initial_value(double q1, double q2) {
  return 1.0 + 0.10 * std::sin(2.0 * pi * q1) * std::cos(2.0 * pi * q2);
}

inline double boundary_value(double q1, double q2) {
  return 0.50 + 0.10 * std::sin(2.0 * pi * q1) * std::cos(2.0 * pi * q2);
}

// Owns all arrays referenced by the non-owning views in Equation and
// BoundaryConditions. Must not be copied or moved after construction.
class CaseData {
 public:
  CaseData(const sayram2d::Mesh& mesh, const CaseConfig& config)
      : mesh_(mesh), config_(config) {
    const std::size_t n1 = mesh_.n1();
    const std::size_t n2 = mesh_.n2();
    const std::size_t n = n1 * n2;
    G_.resize(n);
    D11_.resize(n);
    D12_.resize(n);
    D22_.resize(n);
    inv_tau_.resize(n);
    f_.resize(n);

    for (std::size_t i = 0; i < n1; ++i) {
      const double q1 = mesh_.q1_center(i);
      for (std::size_t j = 0; j < n2; ++j) {
        const double q2 = mesh_.q2_center(j);
        const std::size_t k = i * n2 + j;
        G_[k] = config_.varying_coefficients
                    ? 1.0 + 0.05 * std::sin(2.0 * pi * q1) *
                                std::cos(2.0 * pi * q2)
                    : 1.0;
        D11_[k] = config_.varying_coefficients
                      ? 1.0 + 0.10 * std::cos(2.0 * pi * q1)
                      : 1.0;
        D22_[k] = config_.varying_coefficients
                      ? 0.90 + 0.10 * std::sin(2.0 * pi * q2)
                      : 1.0;
        D12_[k] = config_.off_diagonal ? 0.20 : 0.0;
        inv_tau_[k] = config_.inv_tau;
        f_[k] = initial_value(q1, q2);
      }
    }

    for (std::size_t s = 0; s < 4; ++s) {
      const auto side = static_cast<sayram2d::Side>(s);
      if (config_.bc[s] == sayram2d::BCType::ZeroFlux) {
        bc_.set_zero_flux(side);
        continue;
      }
      if (side == sayram2d::Side::Q1Min ||
          side == sayram2d::Side::Q1Max) {
        const double q1 = (side == sayram2d::Side::Q1Min)
                              ? mesh_.q1_edge(0)
                              : mesh_.q1_edge(mesh_.n1());
        bc_values_[s].resize(n2 + 1);
        for (std::size_t j = 0; j <= n2; ++j) {
          bc_values_[s][j] = boundary_value(q1, mesh_.q2_edge(j));
        }
      } else {
        const double q2 = (side == sayram2d::Side::Q2Min)
                              ? mesh_.q2_edge(0)
                              : mesh_.q2_edge(mesh_.n2());
        bc_values_[s].resize(n1 + 1);
        for (std::size_t i = 0; i <= n1; ++i) {
          bc_values_[s][i] = boundary_value(mesh_.q1_edge(i), q2);
        }
      }
      bc_.set_dirichlet(side, bc_values_[s]);
    }

    equation_ = std::make_unique<sayram2d::Equation>(
        sayram2d::make_view2d(G_, n1, n2),
        sayram2d::make_view2d(D11_, n1, n2),
        sayram2d::make_view2d(D12_, n1, n2),
        sayram2d::make_view2d(D22_, n1, n2), bc_);
    equation_->set_loss(config_.inv_tau);
  }

  CaseData(const CaseData&) = delete;
  CaseData& operator=(const CaseData&) = delete;
  CaseData(CaseData&&) = delete;
  CaseData& operator=(CaseData&&) = delete;

  sayram2d::ArrayView2D<double> f_view() {
    return sayram2d::make_view2d(f_, mesh_.n1(), mesh_.n2());
  }
  sayram2d::ArrayView2D<const double> f_const_view() const {
    return sayram2d::make_view2d(f_, mesh_.n1(), mesh_.n2());
  }
  sayram2d::Equation& equation() { return *equation_; }
  const sayram2d::Equation& equation() const { return *equation_; }
  const std::vector<double>& f() const { return f_; }
  std::vector<double>& f() { return f_; }
  const std::vector<double>& G() const { return G_; }
  std::vector<double>& G() { return G_; }
  std::vector<double>& D11() { return D11_; }
  std::vector<double>& D12() { return D12_; }
  std::vector<double>& D22() { return D22_; }
  const std::vector<double>& inv_tau() const { return inv_tau_; }
  std::vector<double>& inv_tau() { return inv_tau_; }

  void set_boundary_value(sayram2d::Side side, std::size_t index,
                          double value) {
    auto& values = bc_values_[static_cast<std::size_t>(side)];
    if (values.empty()) {
      throw std::invalid_argument(
          "set_boundary_value: side has no array-valued Dirichlet data");
    }
    values.at(index) = value;
  }

  double total_mass(const std::vector<double>& f) const {
    double mass = 0.0;
    for (std::size_t i = 0; i < mesh_.n1(); ++i) {
      for (std::size_t j = 0; j < mesh_.n2(); ++j) {
        const std::size_t k = i * mesh_.n2() + j;
        mass += f[k] * G_[k] * mesh_.cell_measure({i, j});
      }
    }
    return mass;
  }
  double total_mass() const { return total_mass(f_); }

 private:
  const sayram2d::Mesh& mesh_;
  CaseConfig config_;
  std::vector<double> G_;
  std::vector<double> D11_;
  std::vector<double> D12_;
  std::vector<double> D22_;
  std::vector<double> inv_tau_;
  std::vector<double> f_;
  std::array<std::vector<double>, 4> bc_values_;
  sayram2d::BoundaryConditions bc_;
  std::unique_ptr<sayram2d::Equation> equation_;
};

}  // namespace reference
