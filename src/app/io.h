#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "mesh.h"

namespace sayram2d::app {

// Mesh and cell-centered coefficients loaded from the self-describing HDF5
// input file. Matrices are explicitly filled by (i, j).
struct InputData {
  std::size_t n1 = 0;
  std::size_t n2 = 0;
  std::vector<double> q1_edges;
  std::vector<double> q2_edges;
  std::vector<double> q1_centers;
  std::vector<double> q2_centers;
  Eigen::MatrixXd G;
  Eigen::MatrixXd D11;
  Eigen::MatrixXd D12;
  Eigen::MatrixXd D22;
};

InputData load_input(const std::filesystem::path& path);

// Dataset names used by HDF5Output. Defaults are deliberately neutral so the
// writer can be reused by other 2D applications.
struct HDF5OutputNames {
  std::string coordinate1 = "alpha0";
  std::string coordinate2 = "logEN";
  std::string time = "t";
  std::string fields = "PSD";
};

// Native HDF5 writer for a general 2D result layout:
//   /<coordinate1>, /<coordinate2>, /<time>, /<fields>/<frame>
// Matrix datasets are (n1, n2), with the last dimension varying fastest.
class HDF5Output {
 public:
  HDF5Output(const std::filesystem::path& path, std::size_t n1,
             std::size_t n2, HDF5OutputNames names = {});
  ~HDF5Output();

  HDF5Output(const HDF5Output&) = delete;
  HDF5Output& operator=(const HDF5Output&) = delete;

  void write_coordinate1(const std::vector<double>& values);
  void write_coordinate2(const std::vector<double>& values);
  void write_times(const std::vector<double>& values);
  void write_field(std::size_t frame, const Eigen::MatrixXd& field);

 private:
  std::int64_t file_id_ = -1;
  std::int64_t group_id_ = -1;
  std::size_t n1_ = 0;
  std::size_t n2_ = 0;
  HDF5OutputNames names_;
};

}  // namespace sayram2d::app
