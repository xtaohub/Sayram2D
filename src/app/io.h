#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "array2d.h"
#include "mesh.h"

namespace sayram2d::app {

// Mesh and cell-centered coefficients loaded from the self-describing HDF5
// input file. Matrices are explicitly filled by (i, j).
struct InputData {
  std::size_t n1 = 0;
  std::size_t n2 = 0;
  std::vector<double> q1_edges;
  std::vector<double> q2_edges;
  Array2D<double> G;
  Array2D<double> D11;
  Array2D<double> D12;
  Array2D<double> D22;
};

InputData load_input(const std::filesystem::path& path);

// Native HDF5 writer for a general 2D result layout:
//   /<name> for vectors and /<name>/<frame> for fields
// Matrix datasets are (n1, n2), with the last dimension varying fastest.
class HDF5Output {
 public:
  HDF5Output(const std::filesystem::path& path, std::size_t n1,
             std::size_t n2);
  ~HDF5Output();

  HDF5Output(const HDF5Output&) = delete;
  HDF5Output& operator=(const HDF5Output&) = delete;

  // Write a 1D dataset at the root, with an optional `units` attribute.
  void write_vector(const std::vector<double>& values,
                    const std::string& name,
                    const std::string& units = {});

  // Write a 2D field as /<name>/<frame>, with a `units` attribute.
  void write_field(std::size_t frame, ArrayView2D<const double> field,
                   const std::string& name = "PSD",
                   const std::string& units = "arbitrary");

 private:
  std::int64_t file_id_ = -1;
  std::size_t n1_ = 0;
  std::size_t n2_ = 0;
};

}  // namespace sayram2d::app
