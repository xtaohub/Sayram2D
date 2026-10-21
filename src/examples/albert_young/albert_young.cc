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
  const double p = physics::momentum_mev_c(energy);
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

Eigen::MatrixXd AlbertYoungModel::initial_field() const {
  Eigen::MatrixXd field(static_cast<Eigen::Index>(mesh_.n1()),
                        static_cast<Eigen::Index>(mesh_.n2()));
  for (std::size_t i = 0; i < mesh_.n1(); ++i) {
    for (std::size_t j = 0; j < mesh_.n2(); ++j) {
      field(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          albert_young_initial(mesh_.q1_center(i), mesh_.q2_center(j));
    }
  }
  return field;
}

}  // namespace sayram2d::albert_young
