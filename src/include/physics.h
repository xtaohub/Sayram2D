#pragma once

#include <cmath>

namespace sayram2d::physics {

inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double speed_of_light = 299792458.0;  // m/s
inline constexpr double electron_volt = 1.602176634e-19;  // J
inline constexpr double electron_rest_energy_mev = 0.511875;  // MeV
inline constexpr double seconds_per_day = 86400.0;

// Angle conversions.
constexpr double degrees_to_radians(double degrees) noexcept {
  return degrees * pi / 180.0;
}

constexpr double radians_to_degrees(double radians) noexcept {
  return radians * 180.0 / pi;
}

// Energy-unit conversions.
constexpr double ev_to_joule(double ev) noexcept { return ev * electron_volt; }
constexpr double joule_to_ev(double joule) noexcept { return joule / electron_volt; }

constexpr double kev_to_joule(double kev) noexcept { return ev_to_joule(kev * 1.0e3); }
constexpr double joule_to_kev(double joule) noexcept {
  return joule_to_ev(joule) * 1.0e-3;
}

constexpr double mev_to_joule(double mev) noexcept { return ev_to_joule(mev * 1.0e6); }
constexpr double joule_to_mev(double joule) noexcept {
  return joule_to_ev(joule) * 1.0e-6;
}

// Relativistic kinetic-energy/momentum relations.
//
// The three arguments must use a consistent unit system:
//
//   SI:       kinetic_energy and rest_energy in joules,
//             speed_of_light in m/s, momentum in kg m/s;
//   MeV/c:    kinetic_energy and rest_energy in MeV,
//             speed_of_light = 1, momentum in MeV/c.
//
// The caller chooses the unit system; these functions do not assume one.
inline double momentum(double kinetic_energy, double rest_energy,
                       double speed_of_light) noexcept {
  return std::sqrt(kinetic_energy *
                   (kinetic_energy + 2.0 * rest_energy)) /
         speed_of_light;
}

inline double kinetic_energy(double momentum, double rest_energy,
                             double speed_of_light) noexcept {
  const double pc = momentum * speed_of_light;
  return std::sqrt(pc * pc + rest_energy * rest_energy) - rest_energy;
}

}  // namespace sayram2d::physics
