#pragma once

#include <cmath>

#include "physics.h"

namespace sayram2d::radiation_belt {

// Common radiation-belt helper functions for the (alpha0, log E) coordinates.
//
// These helpers are coordinate- and model-specific. They are not used by
// Equation or Solver; callers remain responsible for supplying the G field.
//
// Conventions:
//   alpha0          equatorial pitch angle, radians
//   log_e           log(E / MeV), where E is kinetic energy
//   p               momentum, MeV/c
//   rest_energy_mev particle rest energy, MeV

// d log(E) / dp.
inline double dloge_dp(double log_e, double rest_energy_mev) noexcept {
  const double energy = std::exp(log_e);
  const double p = physics::momentum_mev_c(energy, rest_energy_mev);
  return p / (energy * (energy + rest_energy_mev));
}

// Jacobian for (alpha0, p) coordinates.
inline double G_alpha0_p(double alpha0, double p) noexcept {
  const double sin_alpha0 = std::sin(alpha0);
  return p * p * (1.30 - 0.56 * sin_alpha0) * sin_alpha0 *
         std::cos(alpha0);
}

// Jacobian for (alpha0, log E) coordinates.
inline double G_alpha0_loge(double alpha0, double log_e,
                            double rest_energy_mev) noexcept {
  const double energy = std::exp(log_e);
  const double p = physics::momentum_mev_c(energy, rest_energy_mev);
  return G_alpha0_p(alpha0, p) / dloge_dp(log_e, rest_energy_mev);
}

}  // namespace sayram2d::radiation_belt
