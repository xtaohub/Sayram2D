#include "sayram2d/boundary_conditions.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>

namespace sayram2d {
namespace {

std::size_t expected_size(const Mesh& mesh, Side side) {
  switch (side) {
    case Side::Q1Min:
    case Side::Q1Max:
      return mesh.n2() + 1;
    case Side::Q2Min:
    case Side::Q2Max:
      return mesh.n1() + 1;
  }
  return 0;
}

bool nearly_equal(double a, double b) {
  const double scale = std::max({1.0, std::abs(a), std::abs(b)});
  return std::abs(a - b) <= 1e-12 * scale;
}

}  // namespace

BoundaryConditions::BoundaryConditions() noexcept = default;

void BoundaryConditions::set_dirichlet(Side side, double value) {
  SideData& d = data(side);
  d.type = BCType::Dirichlet;
  d.constant = true;
  d.constant_value = value;
  d.values = {};
}

void BoundaryConditions::set_dirichlet(Side side,
                                       std::span<const double> values) {
  SideData& d = data(side);
  d.type = BCType::Dirichlet;
  d.constant = false;
  d.constant_value = 0.0;
  d.values = values;
}

void BoundaryConditions::set_dirichlet(Side side,
                                       const std::vector<double>& values) {
  set_dirichlet(side, std::span<const double>(values.data(), values.size()));
}

void BoundaryConditions::set_zero_flux(Side side) noexcept {
  SideData& d = data(side);
  d.type = BCType::ZeroFlux;
  d.constant = true;
  d.constant_value = 0.0;
  d.values = {};
}

double BoundaryConditions::constant_value(Side side) const {
  return data(side).constant_value;
}

std::span<const double> BoundaryConditions::values(Side side) const noexcept {
  return data(side).values;
}

double BoundaryConditions::value(Side side, std::size_t index) const {
  const SideData& d = data(side);
  if (d.constant) return d.constant_value;
  assert(index < d.values.size());
  return d.values[index];
}

void BoundaryConditions::validate(const Mesh& mesh) const {
  for (Side side : all_sides()) {
    const SideData& d = data(side);
    if (d.type != BCType::Dirichlet) continue;

    if (d.constant) {
      if (!std::isfinite(d.constant_value)) {
        throw std::invalid_argument("Dirichlet boundary value is not finite");
      }
    } else {
      const std::size_t expected = expected_size(mesh, side);
      if (d.values.size() != expected) {
        throw std::invalid_argument(
            "Dirichlet boundary data has the wrong length for this mesh");
      }
      for (double value : d.values) {
        if (!std::isfinite(value)) {
          throw std::invalid_argument("Dirichlet boundary value is not finite");
        }
      }
    }
  }

  // At a corner, two Dirichlet sides must agree. Dirichlet + ZeroFlux uses the
  // Dirichlet value; ZeroFlux + ZeroFlux keeps the cell-based interpolation.
  for (std::size_t i : {std::size_t{0}, mesh.n1()}) {
    for (std::size_t j : {std::size_t{0}, mesh.n2()}) {
      const Side side1 = (i == 0) ? Side::Q1Min : Side::Q1Max;
      const Side side2 = (j == 0) ? Side::Q2Min : Side::Q2Max;
      const SideData& d1 = data(side1);
      const SideData& d2 = data(side2);
      if (d1.type != BCType::Dirichlet || d2.type != BCType::Dirichlet) {
        continue;
      }
      const double value1 = d1.constant ? d1.constant_value : d1.values[j];
      const double value2 = d2.constant ? d2.constant_value : d2.values[i];
      if (!nearly_equal(value1, value2)) {
        throw std::invalid_argument(
            "inconsistent Dirichlet values at a mesh corner");
      }
    }
  }
}

}  // namespace sayram2d
