// Standalone SAYRAM2D Albert-Young simulation application.
//
// Usage: sayram2d [config.ini]
// The default configuration path is app.ini, resolved from the current
// working directory.

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "albert_young.h"
#include "config.h"
#include "equation.h"
#include "io.h"
#include "mesh.h"
#include "physics.h"
#include "solver.h"

int main(int argc, char** argv) {
  if (argc > 2) {
    std::cerr << "usage: " << argv[0] << " [config.ini]\n";
    return 2;
  }

  const std::filesystem::path config_path =
      (argc == 2) ? argv[1] : "app.ini";

  try {
    const sayram2d::app::Config config =
        sayram2d::app::load_config(config_path);
    const sayram2d::app::InputData input =
        sayram2d::app::load_input(config.data_file);

    sayram2d::Mesh mesh(input.q1_edges, input.q2_edges);
    if (mesh.n1() != input.n1 || mesh.n2() != input.n2) {
      throw std::runtime_error("input dimensions do not match mesh");
    }

    sayram2d::albert_young::AlbertYoungModel model(mesh);
    sayram2d::Equation equation(input.G, input.D11, input.D12, input.D22,
                                model.bc());
    sayram2d::Solver solver(mesh);

    const std::size_t n1 = mesh.n1();
    const std::size_t n2 = mesh.n2();
    const int interval = config.nsteps / (config.n_frames - 1);

    std::vector<double> alpha0_deg(n1);
    std::vector<double> log_en(n2);
    std::vector<double> times(static_cast<std::size_t>(config.n_frames));
    for (std::size_t i = 0; i < n1; ++i) {
      alpha0_deg[i] = sayram2d::physics::radians_to_degrees(mesh.q1_center(i));
    }
    for (std::size_t j = 0; j < n2; ++j) {
      log_en[j] = mesh.q2_center(j) -
                    std::log(sayram2d::physics::electron_rest_energy_mev);
    }
    for (int frame = 0; frame < config.n_frames; ++frame) {
      times[static_cast<std::size_t>(frame)] =
          static_cast<double>(frame * interval) * config.dt;
    }

    const std::filesystem::path output_dir =
        config.output_dir / config.run_id;
    std::filesystem::create_directories(output_dir);
    const std::filesystem::path h5_path =
        output_dir / (config.run_id + "_data.h5");

    const sayram2d::app::HDF5OutputNames output_names{
        "alpha0", "logEN", "t", "f"};
    sayram2d::app::HDF5Output output(h5_path, n1, n2, output_names);
    output.write_coordinate1(alpha0_deg);
    output.write_coordinate2(log_en);
    output.write_times(times);

    Eigen::MatrixXd f_in = model.initial_field();
    Eigen::MatrixXd f_out(static_cast<Eigen::Index>(n1),
                          static_cast<Eigen::Index>(n2));
    output.write_field(0, f_in);

    int saved_frame = 1;
    for (int step = 1; step <= config.nsteps; ++step) {
      solver.solve(equation, f_in, f_out, config.dt);
      std::swap(f_in, f_out);
      if (step % interval == 0) {
        output.write_field(static_cast<std::size_t>(saved_frame), f_in);
        ++saved_frame;
      }
    }

    if (saved_frame != config.n_frames) {
      throw std::runtime_error("internal frame-count mismatch");
    }

    std::cout << "wrote " << h5_path << "\n"
              << "  n1=" << n1 << " n2=" << n2
              << " nsteps=" << config.nsteps << " n_frames="
              << config.n_frames << " dt=" << config.dt << " days\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "sayram2d: " << error.what() << "\n";
    return 1;
  }
}
