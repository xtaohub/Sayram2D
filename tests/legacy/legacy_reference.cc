// Generate reproducible reference outputs for the pre-redesign solver.
//
// Build (macOS, GCC 16, Eigen/xtensor under /Users/xtao/local/include and
// HDF5 headers from Homebrew):
//
//   g++-16 -std=c++17 -O2 \
//     -I/opt/homebrew/include -I/Users/xtao/local/include \
//     -Isource -Isource/Cases \
//     tests/legacy/legacy_reference.cc source/Mesh.cc source/Solver.cc \
//     -o build/legacy_reference
//
// Run:
//
//   ./build/legacy_reference tests/legacy/reference

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Solver.h"
#include "legacy_case.h"

namespace {

using legacy::CaseConfig;
using legacy::LegacyCase;

std::string bc_string(const CaseConfig& config) {
  const char* names[4] = {"Q1Min", "Q1Max", "Q2Min", "Q2Max"};
  std::string result;
  for (std::size_t s = 0; s < 4; ++s) {
    if (!result.empty()) result += ",";
    result += names[s];
    result += ":";
    result += (config.bc[s] == BCType::Dirichlet) ? "Dirichlet" : "ZeroFlux";
  }
  return result;
}

std::vector<double> snapshot(const Solver& solver, const Mesh& mesh) {
  std::vector<double> values(mesh.nx() * mesh.ny());
  for (std::size_t i = 0; i < mesh.nx(); ++i) {
    for (std::size_t j = 0; j < mesh.ny(); ++j) {
      values[j * mesh.nx() + i] = solver.f()(i, j);
    }
  }
  return values;
}

void write_edges(std::ostream& out, const std::vector<double>& edges) {
  out << "[";
  for (std::size_t k = 0; k < edges.size(); ++k) {
    if (k != 0) out << ",";
    out << edges[k];
  }
  out << "]";
}

void run_case(const std::string& name, const std::vector<double>& q1_edges,
              const std::vector<double>& q2_edges, const CaseConfig& config,
              double dt, int steps, const std::string& output_dir) {
  Grid2D grid(q1_edges, q2_edges);
  Mesh mesh(grid, dt);
  LegacyCase equation(mesh, config);
  Solver solver(mesh, &equation);

  std::vector<std::vector<double>> frames;
  frames.reserve(static_cast<std::size_t>(steps) + 1);
  frames.push_back(snapshot(solver, mesh));

  for (int step = 0; step < steps; ++step) {
    solver.update();
    frames.push_back(snapshot(solver, mesh));
  }

  auto mass_of = [&](const std::vector<double>& f) {
    double mass = 0.0;
    for (std::size_t i = 0; i < mesh.nx(); ++i) {
      for (std::size_t j = 0; j < mesh.ny(); ++j) {
        mass += f[j * mesh.nx() + i] * equation.G({i, j}) * mesh.dx(i) *
                mesh.dy(j);
      }
    }
    return mass;
  };
  auto min_of = [&](const std::vector<double>& f) {
    return *std::min_element(f.begin(), f.end());
  };

  const double initial_mass = mass_of(frames.front());
  const double final_mass = mass_of(frames.back());
  const double initial_min = min_of(frames.front());
  const double final_min = min_of(frames.back());

  const std::filesystem::path path =
      std::filesystem::path(output_dir) / (name + ".csv");
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path.string());
  out << std::setprecision(17);

  out << "# case=" << name << "\n";
  out << "# n1=" << mesh.nx() << "\n";
  out << "# n2=" << mesh.ny() << "\n";
  out << "# dt=" << dt << "\n";
  out << "# steps=" << steps << "\n";
  out << "# bc=" << bc_string(config) << "\n";
  out << "# inv_tau=" << config.inv_tau << "\n";
  out << "# off_diagonal=" << (config.off_diagonal ? 1 : 0) << "\n";
  out << "# varying_coefficients="
      << (config.varying_coefficients ? 1 : 0) << "\n";
  out << "# q1_edges=";
  write_edges(out, q1_edges);
  out << "\n# q2_edges=";
  write_edges(out, q2_edges);
  out << "\n";
  out << "# mass0=" << initial_mass << "\n";
  out << "# massN=" << final_mass << "\n";
  out << "# min0=" << initial_min << "\n";
  out << "# minN=" << final_min << "\n";

  out << "i,j,q1,q2,G,D11,D12,D22,inv_tau";
  for (int step = 0; step <= steps; ++step) out << ",f" << step;
  out << "\n";

  for (std::size_t i = 0; i < mesh.nx(); ++i) {
    for (std::size_t j = 0; j < mesh.ny(); ++j) {
      const Ind ind{i, j};
      out << i << "," << j << "," << mesh.x(i) << "," << mesh.y(j) << ","
          << equation.G(ind) << "," << equation.Dxx(ind) << ","
          << equation.Dxy(ind) << "," << equation.Dyy(ind) << ","
          << equation.inv_tau(ind);
      for (int step = 0; step <= steps; ++step) {
        out << "," << frames[static_cast<std::size_t>(step)][j * mesh.nx() + i];
      }
      out << "\n";
    }
  }

  std::cout << name << " n1=" << mesh.nx() << " n2=" << mesh.ny()
            << " mass0=" << initial_mass << " massN=" << final_mass
            << " min0=" << initial_min << " minN=" << final_min << "\n";
}

}  // namespace

int main(int argc, char** argv) {
  const std::string output_dir =
      (argc > 1) ? argv[1] : "tests/legacy/reference";
  std::filesystem::create_directories(output_dir);

  const auto uniform_q1 = legacy::uniform_edges(0.0, 1.0, 12);
  const auto uniform_q2 = legacy::uniform_edges(0.0, 1.2, 10);
  const auto nonuniform_q1 = legacy::nonuniform_q1_edges();
  const auto nonuniform_q2 = legacy::nonuniform_q2_edges();

  const double dt = 0.01;
  const int steps = 5;

  {
    CaseConfig config;
    config.bc = legacy::all_dirichlet();
    config.inv_tau = 0.0;
    config.off_diagonal = false;
    config.varying_coefficients = true;
    run_case("uniform_diag_all_dirichlet", uniform_q1, uniform_q2, config, dt,
             steps, output_dir);
  }

  {
    CaseConfig config;
    config.bc = legacy::mixed_dirichlet_zeroflux();
    config.inv_tau = 0.0;
    config.off_diagonal = false;
    config.varying_coefficients = true;
    run_case("nonuniform_diag_mixed", nonuniform_q1, nonuniform_q2, config, dt,
             steps, output_dir);
  }

  {
    CaseConfig config;
    config.bc = legacy::mixed_dirichlet_zeroflux();
    config.inv_tau = 0.0;
    config.off_diagonal = true;
    config.varying_coefficients = true;
    run_case("nonuniform_offdiag_mixed", nonuniform_q1, nonuniform_q2, config,
             dt, steps, output_dir);
  }

  {
    CaseConfig config;
    config.bc = legacy::all_zeroflux();
    config.inv_tau = 0.0;
    config.off_diagonal = true;
    config.varying_coefficients = true;
    run_case("uniform_offdiag_all_zeroflux", uniform_q1, uniform_q2, config, dt,
             steps, output_dir);
  }

  {
    CaseConfig config;
    config.bc = legacy::mixed_dirichlet_zeroflux();
    config.inv_tau = 0.5;
    config.off_diagonal = true;
    config.varying_coefficients = true;
    run_case("nonuniform_offdiag_loss", nonuniform_q1, nonuniform_q2, config, dt,
             steps, output_dir);
  }

  {
    // The legacy implementation has no explicit source branch. This case
    // records its behavior for negative inv_tau; it is expected to differ
    // from the algorithm document's source treatment.
    CaseConfig config;
    config.bc = legacy::mixed_dirichlet_zeroflux();
    config.inv_tau = -0.1;
    config.off_diagonal = true;
    config.varying_coefficients = true;
    run_case("nonuniform_offdiag_source_legacy", nonuniform_q1, nonuniform_q2,
             config, dt, steps, output_dir);
  }

  return 0;
}
