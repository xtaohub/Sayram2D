#include "io.h"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

namespace sayram2d::app {
namespace {

[[noreturn]] void input_error(const std::filesystem::path& path,
                              const std::string& message) {
  throw std::runtime_error("input " + path.string() + ": " + message);
}

void read_bytes(std::ifstream& stream, void* data, std::size_t bytes,
                const std::filesystem::path& path, const std::string& name) {
  stream.read(static_cast<char*>(data), static_cast<std::streamsize>(bytes));
  if (!stream) input_error(path, "truncated while reading " + name);
}

std::vector<double> read_vector(std::ifstream& stream, std::size_t count,
                                const std::filesystem::path& path,
                                const std::string& name) {
  std::vector<double> values(count);
  read_bytes(stream, values.data(), count * sizeof(double), path, name);
  return values;
}

Eigen::MatrixXd read_matrix(std::ifstream& stream, std::size_t n1,
                            std::size_t n2, const std::filesystem::path& path,
                            const std::string& name) {
  const std::size_t count = n1 * n2;
  std::vector<double> row_major(count);
  read_bytes(stream, row_major.data(), count * sizeof(double), path, name);
  Eigen::MatrixXd matrix(static_cast<Eigen::Index>(n1),
                         static_cast<Eigen::Index>(n2));
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      matrix(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          row_major[i * n2 + j];
    }
  }
  return matrix;
}

void check_status(herr_t status, const std::string& operation) {
  if (status < 0) throw std::runtime_error("HDF5: " + operation + " failed");
}

hid_t check_id(hid_t id, const std::string& operation) {
  if (id < 0) throw std::runtime_error("HDF5: " + operation + " failed");
  return id;
}

class Space {
 public:
  explicit Space(hid_t id) : id_(id) {
    if (id_ < 0) throw std::runtime_error("HDF5: dataspace creation failed");
  }
  ~Space() { H5Sclose(id_); }
  Space(const Space&) = delete;
  Space& operator=(const Space&) = delete;
  hid_t id() const noexcept { return id_; }

 private:
  hid_t id_;
};

class Dataset {
 public:
  explicit Dataset(hid_t id) : id_(id) {
    if (id_ < 0) throw std::runtime_error("HDF5: dataset creation failed");
  }
  ~Dataset() { H5Dclose(id_); }
  Dataset(const Dataset&) = delete;
  Dataset& operator=(const Dataset&) = delete;
  hid_t id() const noexcept { return id_; }

 private:
  hid_t id_;
};

void write_vector_dataset(hid_t file, const char* name,
                          const std::vector<double>& values) {
  const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
  Space space(check_id(H5Screate_simple(1, dims, nullptr), "create dataspace"));
  Dataset dataset(check_id(
      H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space.id(), H5P_DEFAULT,
                 H5P_DEFAULT, H5P_DEFAULT),
      "create dataset " + std::string(name)));
  check_status(H5Dwrite(dataset.id(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, values.data()),
               "write dataset " + std::string(name));
}

}  // namespace

InputData load_input(const std::filesystem::path& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) input_error(path, "cannot open file");

  char magic[8] = {};
  std::int32_t version = 0;
  std::int32_t n1 = 0;
  std::int32_t n2 = 0;
  read_bytes(stream, magic, sizeof(magic), path, "magic");
  read_bytes(stream, &version, sizeof(version), path, "version");
  read_bytes(stream, &n1, sizeof(n1), path, "n1");
  read_bytes(stream, &n2, sizeof(n2), path, "n2");

  constexpr char expected_magic[] = "S2DAYD1";
  if (std::string(magic, 7) != expected_magic || magic[7] != '\0' ||
      version != 1) {
    input_error(path, "invalid S2DAYD1 header");
  }
  if (n1 <= 0 || n2 <= 0) {
    input_error(path, "n1 and n2 must be positive");
  }

  InputData data;
  data.n1 = static_cast<std::size_t>(n1);
  data.n2 = static_cast<std::size_t>(n2);
  if (data.n1 > std::numeric_limits<std::size_t>::max() / data.n2) {
    input_error(path, "n1*n2 overflows size_t");
  }

  data.q1_edges = read_vector(stream, data.n1 + 1, path, "q1_edges");
  data.q2_edges = read_vector(stream, data.n2 + 1, path, "q2_edges");
  data.q1_centers = read_vector(stream, data.n1, path, "q1_centers");
  data.q2_centers = read_vector(stream, data.n2, path, "q2_centers");
  data.G = read_matrix(stream, data.n1, data.n2, path, "G");
  data.D11 = read_matrix(stream, data.n1, data.n2, path, "D11");
  data.D12 = read_matrix(stream, data.n1, data.n2, path, "D12");
  data.D22 = read_matrix(stream, data.n1, data.n2, path, "D22");

  if (stream.peek() != std::char_traits<char>::eof()) {
    input_error(path, "trailing data after coefficient arrays");
  }
  return data;
}

HDF5Output::HDF5Output(const std::filesystem::path& path, std::size_t n1,
                       std::size_t n2)
    : n1_(n1), n2_(n2) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  const hid_t file = check_id(
      H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                H5P_DEFAULT),
      "create " + path.string());
  file_id_ = file;
  try {
    group_id_ = check_id(H5Gcreate2(file_id_, "f", H5P_DEFAULT, H5P_DEFAULT,
                                    H5P_DEFAULT),
                         "create group /f");
  } catch (...) {
    H5Fclose(file_id_);
    file_id_ = -1;
    throw;
  }
}

HDF5Output::~HDF5Output() {
  if (group_id_ >= 0) H5Gclose(group_id_);
  if (file_id_ >= 0) H5Fclose(file_id_);
}

void HDF5Output::write_alpha0(const std::vector<double>& degrees) {
  if (degrees.size() != n1_) {
    throw std::runtime_error("HDF5: /alpha0 has the wrong length");
  }
  write_vector_dataset(file_id_, "alpha0", degrees);
}

void HDF5Output::write_log_en(const std::vector<double>& values) {
  if (values.size() != n2_) {
    throw std::runtime_error("HDF5: /logEN has the wrong length");
  }
  write_vector_dataset(file_id_, "logEN", values);
}

void HDF5Output::write_times(const std::vector<double>& days) {
  write_vector_dataset(file_id_, "t", days);
}

void HDF5Output::write_field(std::size_t frame, const Eigen::MatrixXd& field) {
  if (field.rows() != static_cast<Eigen::Index>(n1_) ||
      field.cols() != static_cast<Eigen::Index>(n2_)) {
    throw std::runtime_error("HDF5: /f field has the wrong dimensions");
  }

  std::vector<double> row_major(n1_ * n2_);
  for (std::size_t i = 0; i < n1_; ++i) {
    for (std::size_t j = 0; j < n2_; ++j) {
      row_major[i * n2_ + j] =
          field(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
    }
  }

  const hsize_t dims[2] = {static_cast<hsize_t>(n1_),
                           static_cast<hsize_t>(n2_)};
  Space space(check_id(H5Screate_simple(2, dims, nullptr), "create dataspace"));
  const std::string name = "f/" + std::to_string(frame);
  Dataset dataset(check_id(
      H5Dcreate2(file_id_, name.c_str(), H5T_NATIVE_DOUBLE, space.id(),
                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
      "create dataset /" + name));
  check_status(H5Dwrite(dataset.id(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, row_major.data()),
               "write dataset /" + name);
  check_status(H5Fflush(file_id_, H5F_SCOPE_GLOBAL), "flush output file");
}

}  // namespace sayram2d::app
