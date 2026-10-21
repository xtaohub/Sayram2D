#pragma once

#include <cmath>

namespace sayram2d::physics {

inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double speed_of_light = 299792458.0;  // m/s
inline constexpr double c = speed_of_light;
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
// The energy arguments are kinetic energies; rest_energy is the particle rest
// energy. SI versions use joules and return kg m/s (or joules).
inline double momentum_si(double kinetic_energy_j, double rest_energy_j) noexcept {
  return std::sqrt(kinetic_energy_j * (kinetic_energy_j + 2.0 * rest_energy_j)) /
         speed_of_light;
}

inline double kinetic_energy_si(double momentum, double rest_energy_j) noexcept {
  const double pc = momentum * speed_of_light;
  return std::sqrt(pc * pc + rest_energy_j * rest_energy_j) - rest_energy_j;
}

// Convenience for callers working directly in electron-volt units.
inline double momentum_si_from_ev(double kinetic_energy_ev,
                                  double rest_energy_ev) noexcept {
  return momentum_si(ev_to_joule(kinetic_energy_ev),
                     ev_to_joule(rest_energy_ev));
}

inline double kinetic_energy_si_from_ev(double momentum,
                                        double rest_energy_ev) noexcept {
  return kinetic_energy_si(momentum, ev_to_joule(rest_energy_ev));
}

// MeV/c momentum for kinetic and rest energies expressed in MeV.
inline double momentum_mev_c(double kinetic_energy_mev,
                             double rest_energy_mev =
                                 electron_rest_energy_mev) noexcept {
  return std::sqrt(kinetic_energy_mev *
                   (kinetic_energy_mev + 2.0 * rest_energy_mev));
}

inline double kinetic_energy_mev_from_momentum_mev_c(
    double momentum_mev_c,
    double rest_energy_mev = electron_rest_energy_mev) noexcept {
  return std::sqrt(momentum_mev_c * momentum_mev_c +
                   rest_energy_mev * rest_energy_mev) -
         rest_energy_mev;
}

}  // namespace sayram2d::physics
