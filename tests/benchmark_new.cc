// Benchmark the redesigned solver. The timing counters are provided by a
// temporary instrumented copy of src/solver.cpp; numerical code is identical
// apart from std::chrono calls.
//
// Usage: benchmark_new <n1> <n2> <steps> [results.csv]

#include <sys/resource.h>

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "reference_case.h"
#include "sayram2d/solver.h"

extern double g_new_vertex_sec;
extern double g_new_assemble_sec;
extern double g_new_factorize_sec;
extern double g_new_solve_sec;
extern double g_new_scatter_sec;

namespace {

using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}

double max_rss_bytes() {
  struct rusage usage {};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return -1.0;
  return static_cast<double>(usage.ru_maxrss);  // macOS: bytes
}

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
      (argc == 5) ? argv[4] : "tests/results/benchmark_new.csv";
  if (n1 == 0 || n2 == 0 || steps <= 0) {
    std::cerr << "n1, n2, and steps must be positive\n";
    return 2;
  }

  const double dt = 0.01;
  reference::CaseConfig config;
  config.bc = reference::mixed_dirichlet_zeroflux();
  config.inv_tau = 0.0;
  config.off_diagonal = false;
  config.varying_coefficients = false;

  const auto q1_edges = reference::uniform_edges(0.0, 1.0, n1);
  const auto q2_edges = reference::uniform_edges(0.0, 1.2, n2);

  const auto construction_start = Clock::now();
  sayram2d::Mesh mesh(q1_edges, q2_edges);
  reference::CaseData data(mesh, config);
  sayram2d::Solver solver(mesh);
  const auto construction_end = Clock::now();
  const double construction = seconds(construction_start, construction_end);

  double first_total = 0.0;
  double first_vertex = 0.0;
  double first_assemble = 0.0;
  double first_factorize = 0.0;
  double first_solve = 0.0;
  double first_scatter = 0.0;
  double repeated_total = 0.0;
  double repeated_vertex = 0.0;
  double repeated_assemble = 0.0;
  double repeated_factorize = 0.0;
  double repeated_solve = 0.0;
  double repeated_scatter = 0.0;

  for (int step = 0; step < steps; ++step) {
    const double v0 = g_new_vertex_sec;
    const double a0 = g_new_assemble_sec;
    const double f0 = g_new_factorize_sec;
    const double s0 = g_new_solve_sec;
    const double c0 = g_new_scatter_sec;

    const auto start = Clock::now();
    solver.solve(data.equation(), data.f_const_view(), data.f_view(), dt);
    const auto end = Clock::now();

    const double total = seconds(start, end);
    const double vertex = g_new_vertex_sec - v0;
    const double assemble = g_new_assemble_sec - a0;
    const double factorize = g_new_factorize_sec - f0;
    const double linear_solve = g_new_solve_sec - s0;
    const double scatter = g_new_scatter_sec - c0;

    if (step == 0) {
      first_total = total;
      first_vertex = vertex;
      first_assemble = assemble;
      first_factorize = factorize;
      first_solve = linear_solve;
      first_scatter = scatter;
    } else {
      repeated_total += total;
      repeated_vertex += vertex;
      repeated_assemble += assemble;
      repeated_factorize += factorize;
      repeated_solve += linear_solve;
      repeated_scatter += scatter;
    }
  }

  if (steps > 1) {
    const double denom = static_cast<double>(steps - 1);
    repeated_total /= denom;
    repeated_vertex /= denom;
    repeated_assemble /= denom;
    repeated_factorize /= denom;
    repeated_solve /= denom;
    repeated_scatter /= denom;
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
           "first_vertex_sec,first_assemble_sec,first_factorize_sec,"
           "first_solve_sec,first_scatter_sec,repeated_total_sec,"
           "repeated_vertex_sec,repeated_assemble_sec,"
           "repeated_factorize_sec,repeated_solve_sec,"
           "repeated_scatter_sec,max_rss_bytes\n";
  }
  out << n1 << "," << n2 << "," << steps << "," << dt << "," << construction
      << "," << first_total << "," << first_vertex << "," << first_assemble
      << "," << first_factorize << "," << first_solve << "," << first_scatter
      << "," << repeated_total << "," << repeated_vertex << ","
      << repeated_assemble << "," << repeated_factorize << ","
      << repeated_solve << "," << repeated_scatter << "," << rss << "\n";

  std::cout << n1 << "x" << n2 << " steps=" << steps
            << " construction=" << construction
            << " first_total=" << first_total
            << " first_factorize=" << first_factorize
            << " first_solve=" << first_solve
            << " repeated_total=" << repeated_total << " rss=" << rss << "\n";
  return 0;
}
