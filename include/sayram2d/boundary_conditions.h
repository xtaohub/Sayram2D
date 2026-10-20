#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "sayram2d/mesh.h"

namespace sayram2d {

enum class BCType { Dirichlet, ZeroFlux };

// Boundary data are non-owning. Array-valued Dirichlet data must remain valid
// until the BoundaryConditions object (and any Equation that copied it) is no
// longer used. Values may be updated in place between solves.
class BoundaryConditions {
 public:
  BoundaryConditions() noexcept;

  void set_dirichlet(Side side, double value);
  void set_dirichlet(Side side, std::span<const double> values);
  void set_dirichlet(Side side, const std::vector<double>& values);
  void set_dirichlet(Side side, std::vector<double>&& values) = delete;

  void set_zero_flux(Side side) noexcept;

  BCType type(Side side) const noexcept { return data(side).type; }
  bool is_dirichlet(Side side) const noexcept {
    return type(side) == BCType::Dirichlet;
  }
  bool is_constant(Side side) const noexcept { return data(side).constant; }
  double constant_value(Side side) const;
  std::span<const double> values(Side side) const noexcept;
  double value(Side side, std::size_t index) const;

  // Checks lengths, finiteness, and corner consistency. Throws
  // std::invalid_argument on failure.
  void validate(const Mesh& mesh) const;

 private:
  struct SideData {
    BCType type = BCType::ZeroFlux;
    bool constant = true;
    double constant_value = 0.0;
    std::span<const double> values{};
  };

  static std::size_t index(Side side) noexcept {
    return static_cast<std::size_t>(side);
  }
  const SideData& data(Side side) const noexcept { return sides_[index(side)]; }
  SideData& data(Side side) noexcept { return sides_[index(side)]; }

  std::array<SideData, 4> sides_{};
};

}  // namespace sayram2d
