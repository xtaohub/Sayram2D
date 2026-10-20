#pragma once

// Baseline test case for the pre-redesign SAYRAM2D implementation.
// This header is intentionally tied to the legacy Equation/Mesh API.
// It is used only by the baseline reference and benchmark executables.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Equation.h"
#include "Mesh.h"

namespace legacy {

constexpr double pi = 3.141592653589793238462643383279502884;

struct CaseConfig {
  // Order follows BoundaryID: XMIN, XMAX, YMIN, YMAX, i.e.
  // q1-min, q1-max, q2-min, q2-max.
  std::array<BCType, 4> bc{BCType::ZeroFlux, BCType::ZeroFlux,
                           BCType::ZeroFlux, BCType::ZeroFlux};
  double inv_tau = 0.0;
  bool off_diagonal = false;
  bool varying_coefficients = true;
};

inline std::array<BCType, 4> all_dirichlet() {
  return {BCType::Dirichlet, BCType::Dirichlet,
          BCType::Dirichlet, BCType::Dirichlet};
}

inline std::array<BCType, 4> mixed_dirichlet_zeroflux() {
  return {BCType::Dirichlet, BCType::ZeroFlux,
          BCType::Dirichlet, BCType::ZeroFlux};
}

inline std::array<BCType, 4> all_zeroflux() {
  return {BCType::ZeroFlux, BCType::ZeroFlux,
          BCType::ZeroFlux, BCType::ZeroFlux};
}

inline std::vector<double> uniform_edges(double a, double b, std::size_t n) {
  if (n == 0) throw std::runtime_error("uniform_edges: n must be positive");
  std::vector<double> e(n + 1);
  const double h = (b - a) / static_cast<double>(n);
  for (std::size_t k = 0; k <= n; ++k) e[k] = a + h * static_cast<double>(k);
  return e;
}

inline std::vector<double> nonuniform_q1_edges() {
  return {0.00, 0.03, 0.08, 0.18, 0.32, 0.50, 0.70, 0.85, 0.95, 0.99, 1.00};
}

inline std::vector<double> nonuniform_q2_edges() {
  return {0.00, 0.05, 0.12, 0.22, 0.35, 0.50, 0.65, 0.78, 0.88, 0.95, 1.00, 1.20};
}

// Smooth, strictly positive initial and boundary data. The same function is
// used on every Dirichlet side so that corner values are consistent.
inline double initial_value(double q1, double q2) {
  return 1.0 + 0.10 * std::sin(2.0 * pi * q1) * std::cos(2.0 * pi * q2);
}

inline double boundary_value(double q1, double q2) {
  return 0.50 + 0.10 * std::sin(2.0 * pi * q1) * std::cos(2.0 * pi * q2);
}

class LegacyCase : public Equation {
 public:
  LegacyCase(const Mesh& mesh, const CaseConfig& config)
      : Equation(mesh), mesh_(mesh), config_(config) {
    fill_coefficients();
  }

  BCType bc_type(BoundaryID side) const override {
    return config_.bc.at(static_cast<std::size_t>(side));
  }

  double init_f(const Ind& ind) const override {
    return initial_value(mesh_.x(ind.i), mesh_.y(ind.j));
  }

  void update(double /*t*/) override {}

  bool dirichlet_vertex_value(BoundaryID side, std::size_t i, std::size_t j,
                              double /*t*/, double* out_value) const override {
    if (bc_type(side) != BCType::Dirichlet) return false;
    *out_value = boundary_value(mesh_.x_edge(i), mesh_.y_edge(j));
    return true;
  }

 private:
  void fill_coefficients() {
    for (std::size_t i = 0; i < mesh_.nx(); ++i) {
      const double q1 = mesh_.x(i);
      for (std::size_t j = 0; j < mesh_.ny(); ++j) {
        const double q2 = mesh_.y(j);
        G_(i, j) = config_.varying_coefficients
                       ? 1.0 + 0.05 * std::sin(2.0 * pi * q1) *
                                   std::cos(2.0 * pi * q2)
                       : 1.0;
        Dxx_(i, j) = config_.varying_coefficients
                         ? 1.0 + 0.10 * std::cos(2.0 * pi * q1)
                         : 1.0;
        Dyy_(i, j) = config_.varying_coefficients
                         ? 0.90 + 0.10 * std::sin(2.0 * pi * q2)
                         : 1.0;
        Dxy_(i, j) = config_.off_diagonal ? 0.20 : 0.0;
        inv_tau_(i, j) = config_.inv_tau;
      }
    }
  }

  const Mesh& mesh_;
  CaseConfig config_;
};

inline double total_mass(const Mesh& mesh, const LegacyCase& eq,
                         const Xtensor2d& f) {
  double mass = 0.0;
  for (std::size_t i = 0; i < mesh.nx(); ++i) {
    for (std::size_t j = 0; j < mesh.ny(); ++j) {
      mass += f(i, j) * eq.G({i, j}) * mesh.dx(i) * mesh.dy(j);
    }
  }
  return mass;
}

inline double minimum_value(const Mesh& mesh, const Xtensor2d& f) {
  double value = f(0, 0);
  for (std::size_t i = 0; i < mesh.nx(); ++i) {
    for (std::size_t j = 0; j < mesh.ny(); ++j) {
      value = std::min(value, f(i, j));
    }
  }
  return value;
}

}  // namespace legacy
