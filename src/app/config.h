#pragma once

#include <filesystem>
#include <string>

namespace sayram2d::app {

// Application configuration. Paths are used as supplied, so relative paths
// resolve against the process working directory.
struct Config {
  std::string run_id;
  std::filesystem::path output_dir = "output";
  std::filesystem::path data_file;
  double dt = 0.0;  // days
  int nsteps = 0;
  int n_frames = 0;  // includes initial and final frames
};

// Reads the strict INI schema described in README.md. Throws
// std::runtime_error with the file and line for malformed or missing values.
Config load_config(const std::filesystem::path& path);

}  // namespace sayram2d::app
