// SAYRAM2D Phase 2 example: the Albert-Young 80x80 case.
//
// The converted data file is produced by tools/convert_albert_young.py and
// already contains cell-centered G, D11, D12, D22. This example performs no
// coefficient interpolation or unit conversion; it only loads the arrays,
// sets up the physical initial/boundary data, advances the solver, updates a
// coefficient in place, and writes a CSV field.
//
// Build with:
//   make example
// or with CMake:
//   cmake -S . -B build-cmake && cmake --build build-cmake
//
// Run from the repository root:
//   ./build/sayram2d/example data/albert_young/ALBERT_YOUNG_80x80.bin albert_young_example.csv

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "sayram2d/array_view2d.h"
#include "sayram2d/boundary_conditions.h"
#include "sayram2d/equation.h"
#include "sayram2d/mesh.h"
#include "sayram2d/solver.h"

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double electron_rest_energy = 0.511875;  // MeV
constexpr double epsilon = std::numeric_limits<double>::epsilon();

struct ConvertedData {
  std::size_t n1 = 0;
  std::size_t n2 = 0;
  std::vector<double> q1_edges;
  std::vector<double> q2_edges;
  std::vector<double> q1_centers;
  std::vector<double> q2_centers;
  std::vector<double> G;
  std::vector<double> D11;
  std::vector<double> D12;
  std::vector<double> D22;
};

ConvertedData load_converted_data(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::runtime_error("cannot open " + path);

  char magic[8] = {};
  std::int32_t version = 0;
  std::int32_t n1 = 0;
  std::int32_t n2 = 0;
  stream.read(magic, sizeof(magic));
  stream.read(reinterpret_cast<char*>(&version), sizeof(version));
  stream.read(reinterpret_cast<char*>(&n1), sizeof(n1));
  stream.read(reinterpret_cast<char*>(&n2), sizeof(n2));
  if (!stream || std::string(magic, 7) != "S2DAYD1" || version != 1 ||
      n1 <= 0 || n2 <= 0) {
    throw std::runtime_error("invalid converted data header: " + path);
  }

  ConvertedData data;
  data.n1 = static_cast<std::size_t>(n1);
  data.n2 = static_cast<std::size_t>(n2);
  const std::size_t cells = data.n1 * data.n2;
  auto read_array = [&](std::vector<double>& values, std::size_t count,
                        const char* name) {
    values.resize(count);
    stream.read(reinterpret_cast<char*>(values.data()),
                static_cast<std::streamsize>(count * sizeof(double)));
    if (!stream) {
      throw std::runtime_error(std::string("truncated converted data: ") +
                               name);
    }
  };

  read_array(data.q1_edges, data.n1 + 1, "q1_edges");
  read_array(data.q2_edges, data.n2 + 1, "q2_edges");
  read_array(data.q1_centers, data.n1, "q1_centers");
  read_array(data.q2_centers, data.n2, "q2_centers");
  read_array(data.G, cells, "G");
  read_array(data.D11, cells, "D11");
  read_array(data.D12, cells, "D12");
  read_array(data.D22, cells, "D22");
  return data;
}

double momentum(double energy_mev) {
  return std::sqrt(energy_mev * (energy_mev + 2.0 * electron_rest_energy));
}

// Albert-Young initial distribution in (alpha0, log E) coordinates.
double initial_value(double alpha, double log_energy) {
  const double energy = std::exp(log_energy);
  const double p = momentum(energy);
  return std::exp(-(energy - 0.2) / 0.1) *
             (std::sin(alpha) - std::sin(5.0 * pi / 180.0)) / (p * p) +
         epsilon;
}

void write_csv(const std::string& path, const sayram2d::Mesh& mesh,
               const std::vector<double>& field) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << std::setprecision(17);
  out << "i,j,alpha0_deg,E_MeV,f\n";
  for (std::size_t i = 0; i < mesh.n1(); ++i) {
    for (std::size_t j = 0; j < mesh.n2(); ++j) {
      out << i << "," << j << "," << mesh.q1_center(i) * 180.0 / pi << ","
          << std::exp(mesh.q2_center(j)) << ","
          << field[i * mesh.n2() + j] << "\n";
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::string data_path =
        (argc > 1) ? argv[1]
                   : "data/albert_young/ALBERT_YOUNG_80x80.bin";
    const std::string output_path =
        (argc > 2) ? argv[2] : "albert_young_example.csv";

    ConvertedData data = load_converted_data(data_path);

    // The former Grid2D class was merged into Mesh in the redesigned API.
    // q1_edges/q2_edges are the full grid specification.
    sayram2d::Mesh mesh(data.q1_edges, data.q2_edges);
    const std::size_t n1 = mesh.n1();
    const std::size_t n2 = mesh.n2();
    const std::size_t cells = n1 * n2;

    // Albert-Young boundary conditions:
    //   alpha0 = 5 deg: f = 0 (Dirichlet)
    //   alpha0 = 90 deg: zero flux
    //   E = 0.2 MeV: initial spectrum (Dirichlet along the boundary)
    //   E = 5 MeV: f = 0 (Dirichlet)
    std::vector<double> lower_energy_boundary(n1 + 1);
    for (std::size_t i = 0; i <= n1; ++i) {
      lower_energy_boundary[i] =
          initial_value(mesh.q1_edge(i), mesh.q2_edge(0));
    }

    sayram2d::BoundaryConditions bc;
    bc.set_dirichlet(sayram2d::Side::Q1Min, 0.0);
    bc.set_zero_flux(sayram2d::Side::Q1Max);
    bc.set_dirichlet(sayram2d::Side::Q2Min, lower_energy_boundary);
    bc.set_dirichlet(sayram2d::Side::Q2Max, 0.0);

    // G, D11, D12, D22 are already cell-centered and in simulation units.
    sayram2d::Equation equation(
        sayram2d::make_view2d(data.G, n1, n2),
        sayram2d::make_view2d(data.D11, n1, n2),
        sayram2d::make_view2d(data.D12, n1, n2),
        sayram2d::make_view2d(data.D22, n1, n2), bc);

    std::vector<double> f_in(cells);
    std::vector<double> f_out(cells, 0.0);
    for (std::size_t i = 0; i < n1; ++i) {
      for (std::size_t j = 0; j < n2; ++j) {
        f_in[i * n2 + j] =
            initial_value(mesh.q1_center(i), mesh.q2_center(j));
      }
    }

    sayram2d::Solver solver(mesh);
    const double dt = 0.002;  // days
    const int steps = 10;

    for (int step = 0; step < steps; ++step) {
      if (step == 5) {
        // Demonstrate a coefficient update without rebuilding the solver.
        // Scaling both diagonal components preserves positive definiteness.
        for (double& value : data.D11) value *= 1.01;
        for (double& value : data.D22) value *= 1.01;
        std::cout << "step " << step << ": updated D11 and D22 in place\n";
      }

      solver.solve(equation, sayram2d::make_view2d(f_in, n1, n2),
                   sayram2d::make_view2d(f_out, n1, n2), dt);
      std::swap(f_in, f_out);
    }

    write_csv(output_path, mesh, f_in);

    const auto minmax = std::minmax_element(f_in.begin(), f_in.end());
    std::cout << "wrote " << output_path << "\n"
              << "  cells=" << cells << " steps=" << steps << " dt=" << dt
              << " days\n"
              << "  f range=[" << *minmax.first << ", " << *minmax.second
              << "]\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Albert-Young example failed: " << error.what() << "\n";
    return 1;
  }
}
