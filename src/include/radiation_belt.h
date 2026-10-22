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
//   log_e           log(E), where E is kinetic energy
//   p               momentum in units consistent with E
//   rest_energy    particle rest energy, in the same energy unit as E
//
// The formulas require E, rest_energy, and speed_of_light to use a consistent
// unit system. For example, with E in MeV and p in MeV/c, pass c = 1; with SI
// units, pass energy in joules and c in m/s.

// d log(E) / dp = p c^2 / [E (E + E0)].
inline double dloge_dp(double log_e, double rest_energy,
                       double speed_of_light) noexcept {
  const double energy = std::exp(log_e);
  const double p = physics::momentum(energy, rest_energy, speed_of_light);
  return speed_of_light * speed_of_light * p /
         (energy * (energy + rest_energy));
}

// Jacobian for (alpha0, p) coordinates.
inline double G_alpha0_p(double alpha0, double p) noexcept {
  const double sin_alpha0 = std::sin(alpha0);
  return p * p * (1.30 - 0.56 * sin_alpha0) * sin_alpha0 *
         std::cos(alpha0);
}

// Jacobian for (alpha0, log E) coordinates.
inline double G_alpha0_loge(double alpha0, double log_e,
                            double rest_energy,
                            double speed_of_light) noexcept {
  const double energy = std::exp(log_e);
  const double p = physics::momentum(energy, rest_energy, speed_of_light);
  return G_alpha0_p(alpha0, p) /
         dloge_dp(log_e, rest_energy, speed_of_light);
}

}  // namespace sayram2d::radiation_belt
