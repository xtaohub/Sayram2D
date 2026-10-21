#pragma once

#include <vector>

#include <Eigen/Dense>

#include "boundary_conditions.h"
#include "mesh.h"

namespace sayram2d::albert_young {

// Original Albert-Young initial distribution in (alpha0, log E) coordinates.
double albert_young_initial(double alpha, double log_energy) noexcept;

// Physical setup of the legacy Albert-Young case. The model owns the
// Q2Min boundary values because BoundaryConditions stores a non-owning span.
class AlbertYoungModel {
 public:
  explicit AlbertYoungModel(const Mesh& mesh);

  AlbertYoungModel(const AlbertYoungModel&) = delete;
  AlbertYoungModel& operator=(const AlbertYoungModel&) = delete;

  const BoundaryConditions& bc() const noexcept { return bc_; }
  const std::vector<double>& q2_min_values() const noexcept {
    return q2_min_values_;
  }
  Eigen::MatrixXd initial_field() const;

 private:
  const Mesh& mesh_;
  std::vector<double> q2_min_values_;
  BoundaryConditions bc_;
};

}  // namespace sayram2d::albert_young
