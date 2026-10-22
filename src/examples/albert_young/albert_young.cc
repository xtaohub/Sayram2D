#include "albert_young.h"

#include <cmath>
#include <cstddef>
#include <limits>

#include "physics.h"

namespace sayram2d::albert_young {
namespace {

constexpr double epsilon = std::numeric_limits<double>::epsilon();

}  // namespace

double albert_young_initial(double alpha, double log_energy) noexcept {
  const double energy = std::exp(log_energy);
  // Natural units for this case: E in MeV and p in MeV/c, so c = 1.
  const double p = physics::momentum(
      energy, physics::electron_rest_energy_mev, 1.0);
  return std::exp(-(energy - 0.2) / 0.1) *
             (std::sin(alpha) -
              std::sin(physics::degrees_to_radians(5.0))) /
             (p * p) +
         epsilon;
}

AlbertYoungModel::AlbertYoungModel(const Mesh& mesh)
    : mesh_(mesh), q2_min_values_(mesh.n1() + 1) {
  for (std::size_t i = 0; i <= mesh_.n1(); ++i) {
    q2_min_values_[i] =
        albert_young_initial(mesh_.q1_edge(i), mesh_.q2_edge(0));
  }

  bc_.set_dirichlet(Side::Q1Min, 0.0);
  bc_.set_zero_flux(Side::Q1Max);
  bc_.set_dirichlet(Side::Q2Min, q2_min_values_);
  bc_.set_dirichlet(Side::Q2Max, 0.0);
}

Array2D<double> AlbertYoungModel::initial_field() const {
  Array2D<double> field(mesh_.n1(), mesh_.n2());
  for (std::size_t j = 0; j < mesh_.n2(); ++j) {
    for (std::size_t i = 0; i < mesh_.n1(); ++i) {
      field(i, j) =
          albert_young_initial(mesh_.q1_center(i), mesh_.q2_center(j));
    }
  }
  return field;
}

}  // namespace sayram2d::albert_young
