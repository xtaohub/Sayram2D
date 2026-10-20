// Benchmark the pre-redesign solver.
//
// The timing counters are provided by a temporary instrumented copy of
// source/Solver.cc; numerical code is identical apart from std::chrono calls.
// See tests/legacy/README.md for the build procedure.
//
// Usage:
//   legacy_benchmark <n1> <n2> <steps> [results.csv]
//
// The solver is advanced `steps` times. The first step is reported separately
// because it includes the first numeric factorization; later steps are
// averaged.

#include <sys/resource.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "Solver.h"
#include "legacy_case.h"

extern double g_legacy_assemble_sec;
extern double g_legacy_factorize_sec;
extern double g_legacy_solve_sec;
extern double g_legacy_scatter_sec;
extern double g_legacy_post_sec;

namespace {

using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}

double max_rss_bytes() {
  struct rusage usage {};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return -1.0;
  // macOS reports bytes; Linux reports kilobytes. This benchmark is run on
  // macOS, so bytes are used and documented in the report.
  return static_cast<double>(usage.ru_maxrss);
}

struct Timings {
  double construction = 0.0;
  double first_total = 0.0;
  double first_assemble = 0.0;
  double first_factorize = 0.0;
  double first_solve = 0.0;
  double first_scatter = 0.0;
  double first_post = 0.0;
  double repeated_total = 0.0;
  double repeated_assemble = 0.0;
  double repeated_factorize = 0.0;
  double repeated_solve = 0.0;
  double repeated_scatter = 0.0;
  double repeated_post = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
  if (argc < 4 || argc > 5) {
    std::cerr << "usage: " << argv[0]
              << " <n1> <n2> <steps> [results.csv]\n";
    return 2;
  }
  const std::size_t n1 = static_cast<std::size_t>(std::stoul(argv[1]));
  const std::size_t n2 = static_cast<std::size_t>(std::stoul(argv[2]));
  const int steps = std::stoi(argv[3]);
  const std::string output_path =
      (argc == 5) ? argv[4] : "tests/legacy/results/benchmark_original.csv";

  if (n1 == 0 || n2 == 0 || steps <= 0) {
    std::cerr << "n1, n2, and steps must be positive\n";
    return 2;
  }

  const double dt = 0.01;
  legacy::CaseConfig config;
  config.bc = legacy::mixed_dirichlet_zeroflux();
  config.inv_tau = 0.0;
  config.off_diagonal = false;
  config.varying_coefficients = false;

  const auto q1_edges = legacy::uniform_edges(0.0, 1.0, n1);
  const auto q2_edges = legacy::uniform_edges(0.0, 1.2, n2);
  Grid2D grid(q1_edges, q2_edges);
  Mesh mesh(grid, dt);

  const auto construction_start = Clock::now();
  legacy::LegacyCase equation(mesh, config);
  Solver solver(mesh, &equation);
  const auto construction_end = Clock::now();

  Timings timings;
  timings.construction = seconds(construction_start, construction_end);

  for (int step = 0; step < steps; ++step) {
    const double assemble0 = g_legacy_assemble_sec;
    const double factorize0 = g_legacy_factorize_sec;
    const double solve0 = g_legacy_solve_sec;
    const double scatter0 = g_legacy_scatter_sec;
    const double post0 = g_legacy_post_sec;

    const auto start = Clock::now();
    solver.update();
    const auto end = Clock::now();

    const double total = seconds(start, end);
    const double assemble = g_legacy_assemble_sec - assemble0;
    const double factorize = g_legacy_factorize_sec - factorize0;
    const double linear_solve = g_legacy_solve_sec - solve0;
    const double scatter = g_legacy_scatter_sec - scatter0;
    const double post = g_legacy_post_sec - post0;

    if (step == 0) {
      timings.first_total = total;
      timings.first_assemble = assemble;
      timings.first_factorize = factorize;
      timings.first_solve = linear_solve;
      timings.first_scatter = scatter;
      timings.first_post = post;
    } else {
      timings.repeated_total += total;
      timings.repeated_assemble += assemble;
      timings.repeated_factorize += factorize;
      timings.repeated_solve += linear_solve;
      timings.repeated_scatter += scatter;
      timings.repeated_post += post;
    }
  }

  if (steps > 1) {
    const double denom = static_cast<double>(steps - 1);
    timings.repeated_total /= denom;
    timings.repeated_assemble /= denom;
    timings.repeated_factorize /= denom;
    timings.repeated_solve /= denom;
    timings.repeated_scatter /= denom;
    timings.repeated_post /= denom;
  }

  const double rss = max_rss_bytes();
  const std::filesystem::path path(output_path);
  std::filesystem::create_directories(path.parent_path());
  const bool write_header = !std::filesystem::exists(path);
  std::ofstream out(path, std::ios::app);
  if (!out) {
    std::cerr << "cannot write " << path << "\n";
    return 1;
  }
  out << std::setprecision(10);
  if (write_header) {
    out << "n1,n2,steps,dt,construction_sec,first_total_sec,"
           "first_assemble_sec,first_factorize_sec,first_solve_sec,"
           "first_scatter_sec,first_post_sec,repeated_total_sec,"
           "repeated_assemble_sec,repeated_factorize_sec,"
           "repeated_solve_sec,repeated_scatter_sec,repeated_post_sec,"
           "max_rss_bytes\n";
  }
  out << n1 << "," << n2 << "," << steps << "," << dt << ","
      << timings.construction << "," << timings.first_total << ","
      << timings.first_assemble << "," << timings.first_factorize << ","
      << timings.first_solve << "," << timings.first_scatter << ","
      << timings.first_post << "," << timings.repeated_total << ","
      << timings.repeated_assemble << "," << timings.repeated_factorize << ","
      << timings.repeated_solve << "," << timings.repeated_scatter << ","
      << timings.repeated_post << "," << rss << "\n";

  std::cout << n1 << "x" << n2 << " steps=" << steps
            << " construction=" << timings.construction
            << " first_total=" << timings.first_total
            << " first_factorize=" << timings.first_factorize
            << " first_solve=" << timings.first_solve
            << " repeated_total=" << timings.repeated_total
            << " rss=" << rss << "\n";
  return 0;
}
