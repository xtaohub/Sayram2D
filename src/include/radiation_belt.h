#pragma once

#include <cmath>

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
// The formulas only require E and rest_energy to use the same energy unit.
// For example, log(E/MeV) requires rest_energy in MeV; log(E/GeV) requires
// rest_energy in GeV.
//
// Natural units with c = 1 are used. With E in MeV and rest_energy in MeV,
// p is in MeV/c. In SI units, p = sqrt(E(E + 2 E0)) / c, and dloge_dp gains
// a factor of c.

// d log(E) / dp.
inline double dloge_dp(double log_e, double rest_energy) noexcept {
  const double energy = std::exp(log_e);
  const double p = std::sqrt(energy * (energy + 2.0 * rest_energy));
  return p / (energy * (energy + rest_energy));
}

// Jacobian for (alpha0, p) coordinates.
inline double G_alpha0_p(double alpha0, double p) noexcept {
  const double sin_alpha0 = std::sin(alpha0);
  return p * p * (1.30 - 0.56 * sin_alpha0) * sin_alpha0 *
         std::cos(alpha0);
}

// Jacobian for (alpha0, log E) coordinates.
inline double G_alpha0_loge(double alpha0, double log_e,
                            double rest_energy) noexcept {
  const double energy = std::exp(log_e);
  const double p = std::sqrt(energy * (energy + 2.0 * rest_energy));
  return G_alpha0_p(alpha0, p) / dloge_dp(log_e, rest_energy);
}

}  // namespace sayram2d::radiation_belt
