// Generate reference outputs with the new library for comparison against the
// legacy baseline and for use by the regression tests.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "reference_case.h"
#include "sayram2d/solver.h"

namespace {

using reference::CaseConfig;
using reference::CaseData;

std::string bc_string(const CaseConfig& config) {
  const char* names[4] = {"Q1Min", "Q1Max", "Q2Min", "Q2Max"};
  std::string result;
  for (std::size_t s = 0; s < 4; ++s) {
    if (!result.empty()) result += ",";
    result += names[s];
    result += ":";
    result += (config.bc[s] == sayram2d::BCType::Dirichlet) ? "Dirichlet"
                                                           : "ZeroFlux";
  }
  return result;
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
  sayram2d::Mesh mesh(q1_edges, q2_edges);
  CaseData data(mesh, config);
  sayram2d::Solver solver(mesh);

  std::vector<std::vector<double>> frames;
  frames.reserve(static_cast<std::size_t>(steps) + 1);
  frames.push_back(data.f());

  for (int step = 0; step < steps; ++step) {
    solver.solve(data.equation(), data.f_const_view(), data.f_view(), dt);
    frames.push_back(data.f());
  }

  const double mass0 = data.total_mass(frames.front());
  const double massN = data.total_mass(frames.back());
  const double min0 = *std::min_element(frames.front().begin(),
                                        frames.front().end());
  const double minN =
      *std::min_element(frames.back().begin(), frames.back().end());

  const std::filesystem::path path =
      std::filesystem::path(output_dir) / (name + ".csv");
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path.string());
  out << std::setprecision(17);
  out << "# case=" << name << "\n";
  out << "# n1=" << mesh.n1() << "\n";
  out << "# n2=" << mesh.n2() << "\n";
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
  out << "# mass0=" << mass0 << "\n";
  out << "# massN=" << massN << "\n";
  out << "# min0=" << min0 << "\n";
  out << "# minN=" << minN << "\n";

  out << "i,j,q1,q2,G,D11,D12,D22,inv_tau";
  for (int step = 0; step <= steps; ++step) out << ",f" << step;
  out << "\n";

  for (std::size_t i = 0; i < mesh.n1(); ++i) {
    for (std::size_t j = 0; j < mesh.n2(); ++j) {
      const std::size_t k = i * mesh.n2() + j;
      out << i << "," << j << "," << mesh.q1_center(i) << ","
          << mesh.q2_center(j) << "," << data.G()[k] << ","
          << data.equation().D11(i, j) << ","
          << data.equation().D12(i, j) << ","
          << data.equation().D22(i, j) << ","
          << data.equation().inv_tau(i, j);
      for (int step = 0; step <= steps; ++step) {
        out << "," << frames[static_cast<std::size_t>(step)][k];
      }
      out << "\n";
    }
  }

  std::cout << name << " n1=" << mesh.n1() << " n2=" << mesh.n2()
            << " mass0=" << mass0 << " massN=" << massN << " min0=" << min0
            << " minN=" << minN << "\n";
}

}  // namespace

int main(int argc, char** argv) {
  const std::string output_dir =
      (argc > 1) ? argv[1] : "build/new_reference";
  std::filesystem::create_directories(output_dir);

  const auto uniform_q1 = reference::uniform_edges(0.0, 1.0, 12);
  const auto uniform_q2 = reference::uniform_edges(0.0, 1.2, 10);
  const auto nonuniform_q1 = reference::nonuniform_q1_edges();
  const auto nonuniform_q2 = reference::nonuniform_q2_edges();
  const double dt = 0.01;
  const int steps = 5;

  {
    CaseConfig config;
    config.bc = reference::all_dirichlet();
    config.inv_tau = 0.0;
    config.off_diagonal = false;
    config.varying_coefficients = true;
    run_case("uniform_diag_all_dirichlet", uniform_q1, uniform_q2, config, dt,
             steps, output_dir);
  }
  {
    CaseConfig config;
    config.bc = reference::mixed_dirichlet_zeroflux();
    config.inv_tau = 0.0;
    config.off_diagonal = false;
    config.varying_coefficients = true;
    run_case("nonuniform_diag_mixed", nonuniform_q1, nonuniform_q2, config, dt,
             steps, output_dir);
  }
  {
    CaseConfig config;
    config.bc = reference::mixed_dirichlet_zeroflux();
    config.inv_tau = 0.0;
    config.off_diagonal = true;
    config.varying_coefficients = true;
    run_case("nonuniform_offdiag_mixed", nonuniform_q1, nonuniform_q2, config,
             dt, steps, output_dir);
  }
  {
    CaseConfig config;
    config.bc = reference::all_zeroflux();
    config.inv_tau = 0.0;
    config.off_diagonal = true;
    config.varying_coefficients = true;
    run_case("uniform_offdiag_all_zeroflux", uniform_q1, uniform_q2, config, dt,
             steps, output_dir);
  }
  {
    CaseConfig config;
    config.bc = reference::mixed_dirichlet_zeroflux();
    config.inv_tau = 0.5;
    config.off_diagonal = true;
    config.varying_coefficients = true;
    run_case("nonuniform_offdiag_loss", nonuniform_q1, nonuniform_q2, config, dt,
             steps, output_dir);
  }
  {
    CaseConfig config;
    config.bc = reference::mixed_dirichlet_zeroflux();
    config.inv_tau = -0.1;
    config.off_diagonal = true;
    config.varying_coefficients = true;
    run_case("nonuniform_offdiag_source", nonuniform_q1, nonuniform_q2, config,
             dt, steps, output_dir);
  }

  return 0;
}
