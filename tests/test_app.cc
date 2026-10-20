#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "albert_young.h"
#include "config.h"
#include "io.h"

namespace {

using namespace sayram2d;
using namespace sayram2d::app;

void check(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

bool near(double a, double b, double tolerance = 1e-12) {
  const double scale = std::max({1.0, std::abs(a), std::abs(b)});
  return std::abs(a - b) <= tolerance * scale;
}

template <typename Function>
void expect_throw(Function&& function, const std::string& message) {
  try {
    function();
  } catch (const std::exception&) {
    return;
  }
  throw std::runtime_error(message);
}

void write_text(const std::filesystem::path& path, const std::string& text) {
  std::ofstream output(path);
  check(static_cast<bool>(output), "cannot write test config");
  output << text;
}

void test_config() {
  const auto config = load_config("tests/app_test.ini");
  check(config.run_id == "AppTest", "run_id");
  check(config.output_dir == "build/app_test_output", "output_dir");
  check(config.data_file == "data/albert_young/ALBERT_YOUNG_80x80.h5",
        "data_file");
  check(near(config.dt, 0.002) && config.nsteps == 10 && config.n_frames == 11,
        "time/output settings");

  const auto temporary =
      std::filesystem::temp_directory_path() / "sayram2d_app_config_test.ini";
  write_text(temporary,
             "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
             "[time]\ndt=0.1\nnsteps=4\n");
  const auto defaults = load_config(temporary);
  check(defaults.output_dir == "output", "default output_dir");
  check(defaults.n_frames == 5, "default n_frames");
  std::filesystem::remove(temporary);

  write_text(temporary, "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
                        "[time]\ndt=0.1\nnsteps=4\n[output]\nn_frames=4\n");
  expect_throw([&] { (void)load_config(temporary); },
               "nondivisible n_frames must throw");

  write_text(temporary, "[run]\nrun_id=Test\n[run]\nrun_id=Other\n"
                        "[input]\ndata_file=data.bin\n[time]\ndt=0.1\n"
                        "nsteps=4\n");
  expect_throw([&] { (void)load_config(temporary); },
               "duplicate key must throw");

  write_text(temporary, "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
                        "[time]\ndt=0.1\nnsteps=4\n[extra]\nkey=1\n");
  expect_throw([&] { (void)load_config(temporary); },
               "unknown section must throw");

  write_text(temporary, "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
                        "[time]\ndt=abc\nnsteps=4\n");
  expect_throw([&] { (void)load_config(temporary); },
               "malformed dt must throw");

  write_text(temporary, "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
                        "[time]\ndt=0\nnsteps=4\n");
  expect_throw([&] { (void)load_config(temporary); },
               "nonpositive dt must throw");

  write_text(temporary, "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
                        "[time]\ndt=0.1\nnsteps=0\n");
  expect_throw([&] { (void)load_config(temporary); },
               "nonpositive nsteps must throw");

  write_text(temporary, "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
                        "[time]\ndt=0.1\nnsteps=4\n[output]\nn_frames=1\n");
  expect_throw([&] { (void)load_config(temporary); },
               "n_frames < 2 must throw");

  write_text(temporary, "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
                        "[time]\ndt=0.1\nnsteps=4\nbogus=1\n");
  expect_throw([&] { (void)load_config(temporary); },
               "unknown key must throw");

  write_text(temporary, "[run]\nrun_id=Test\n[input]\ndata_file=data.bin\n"
                        "[time]\ndt=0.1\n");
  expect_throw([&] { (void)load_config(temporary); },
               "missing nsteps must throw");
  std::filesystem::remove(temporary);
  expect_throw([] { (void)load_config("missing-app-config.ini"); },
               "missing config must throw");
}

void test_input_loader() {
  const auto data = load_input("data/albert_young/ALBERT_YOUNG_80x80.h5");
  check(data.n1 == 80 && data.n2 == 80, "input dimensions");
  check(data.q1_edges.size() == 81 && data.q2_edges.size() == 81,
        "input edge lengths");
  check(data.G.rows() == 80 && data.G.cols() == 80 &&
            data.D22.rows() == 80 && data.D22.cols() == 80,
        "input matrix dimensions");
  check(data.G(0, 0) == data.G(0, 0), "input coefficients are finite");

  const auto temporary =
      std::filesystem::temp_directory_path() / "sayram2d_app_input_test.h5";
  {
    HDF5Output incomplete(temporary, data.n1, data.n2);
    incomplete.write_times({0.0});
  }
  expect_throw(
      [&] { (void)load_input(temporary); },
      "HDF5 input without required datasets must throw");
  std::filesystem::remove(temporary);

  write_text(temporary, "not an HDF5 file");
  expect_throw([&] { (void)load_input(temporary); },
               "corrupt HDF5 input must throw");
  std::filesystem::remove(temporary);
  expect_throw([] { (void)load_input("missing-input.h5"); },
               "missing input must throw");

}

void test_albert_young_setup() {
  const auto data = load_input("data/albert_young/ALBERT_YOUNG_80x80.h5");
  const Mesh mesh(data.q1_edges, data.q2_edges);
  const AlbertYoungModel model(mesh);

  check(model.bc().is_dirichlet(Side::Q1Min) &&
            model.bc().constant_value(Side::Q1Min) == 0.0,
        "Q1Min boundary");
  check(model.bc().type(Side::Q1Max) == BCType::ZeroFlux, "Q1Max boundary");
  check(model.bc().is_dirichlet(Side::Q2Max) &&
            model.bc().constant_value(Side::Q2Max) == 0.0,
        "Q2Max boundary");
  check(model.bc().is_dirichlet(Side::Q2Min), "Q2Min boundary type");
  check(model.q2_min_values().size() == mesh.n1() + 1,
        "Q2Min boundary length");

  for (std::size_t i = 0; i <= mesh.n1(); ++i) {
    check(near(model.q2_min_values()[i],
               albert_young_initial(mesh.q1_edge(i), mesh.q2_edge(0))),
          "Q2Min boundary function");
  }
  check(near(model.q2_min_values().front(), 0.0),
        "Q2Min corner value follows the legacy prescription");
  model.bc().validate(mesh);

  const Eigen::MatrixXd initial = model.initial_field();
  check(initial.rows() == 80 && initial.cols() == 80, "initial field shape");
  check(near(initial(0, 0),
             albert_young_initial(mesh.q1_center(0), mesh.q2_center(0))),
        "initial field formula");

  expect_throw(
      [] {
        HDF5Output output("/nonexistent-directory/result.h5", 2, 2);
        (void)output;
      },
      "invalid HDF5 path must throw");
}

}  // namespace

int main() {
  try {
    test_config();
    test_input_loader();
    test_albert_young_setup();
    std::cout << "All standalone application tests passed.\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "app test failed: " << error.what() << "\n";
    return 1;
  }
}
