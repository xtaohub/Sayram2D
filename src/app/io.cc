#include "io.h"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <hdf5.h>

namespace sayram2d::app {
namespace {

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

class Datatype {
 public:
  explicit Datatype(hid_t id) : id_(id) {
    if (id_ < 0) throw std::runtime_error("HDF5: datatype creation failed");
  }
  ~Datatype() { H5Tclose(id_); }
  Datatype(const Datatype&) = delete;
  Datatype& operator=(const Datatype&) = delete;
  hid_t id() const noexcept { return id_; }

 private:
  hid_t id_;
};

class Attribute {
 public:
  explicit Attribute(hid_t id) : id_(id) {
    if (id_ < 0) throw std::runtime_error("HDF5: attribute creation failed");
  }
  ~Attribute() { H5Aclose(id_); }
  Attribute(const Attribute&) = delete;
  Attribute& operator=(const Attribute&) = delete;
  hid_t id() const noexcept { return id_; }

 private:
  hid_t id_;
};

class File {
 public:
  explicit File(hid_t id) : id_(id) {
    if (id_ < 0) throw std::runtime_error("HDF5: cannot open input file");
  }
  ~File() { H5Fclose(id_); }
  File(const File&) = delete;
  File& operator=(const File&) = delete;
  hid_t id() const noexcept { return id_; }

 private:
  hid_t id_;
};

std::vector<double> read_vector(hid_t file, const char* name) {
  Dataset dataset(check_id(H5Dopen2(file, name, H5P_DEFAULT),
                           "open dataset /" + std::string(name)));
  Space space(check_id(H5Dget_space(dataset.id()),
                       "get dataspace for /" + std::string(name)));
  if (H5Sget_simple_extent_ndims(space.id()) != 1) {
    throw std::runtime_error("HDF5: dataset /" + std::string(name) +
                             " must be one-dimensional");
  }
  hsize_t dimensions[1] = {};
  check_status(H5Sget_simple_extent_dims(space.id(), dimensions, nullptr),
               "get dimensions for /" + std::string(name));
  std::vector<double> values(static_cast<std::size_t>(dimensions[0]));
  check_status(H5Dread(dataset.id(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, values.data()),
               "read dataset /" + std::string(name));
  return values;
}

Eigen::MatrixXd read_matrix(hid_t file, const char* name, std::size_t n1,
                            std::size_t n2) {
  Dataset dataset(check_id(H5Dopen2(file, name, H5P_DEFAULT),
                           "open dataset /" + std::string(name)));
  Space space(check_id(H5Dget_space(dataset.id()),
                       "get dataspace for /" + std::string(name)));
  if (H5Sget_simple_extent_ndims(space.id()) != 2) {
    throw std::runtime_error("HDF5: dataset /" + std::string(name) +
                             " must be two-dimensional");
  }
  hsize_t dimensions[2] = {};
  check_status(H5Sget_simple_extent_dims(space.id(), dimensions, nullptr),
               "get dimensions for /" + std::string(name));
  if (dimensions[0] != n1 || dimensions[1] != n2) {
    throw std::runtime_error("HDF5: dataset /" + std::string(name) +
                             " has the wrong shape");
  }

  std::vector<double> row_major(n1 * n2);
  check_status(H5Dread(dataset.id(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                       H5P_DEFAULT, row_major.data()),
               "read dataset /" + std::string(name));
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

void write_string_attribute(hid_t object, const char* name,
                            const std::string& value) {
  Space space(check_id(H5Screate(H5S_SCALAR), "create attribute dataspace"));
  Datatype type(check_id(H5Tcopy(H5T_C_S1), "copy string datatype"));
  check_status(H5Tset_size(type.id(), value.size() + 1),
               "set string attribute size");
  Attribute attribute(check_id(
      H5Acreate2(object, name, type.id(), space.id(), H5P_DEFAULT,
                 H5P_DEFAULT),
      "create attribute " + std::string(name)));
  check_status(H5Awrite(attribute.id(), type.id(), value.c_str()),
               "write attribute " + std::string(name));
}

void write_vector_dataset(hid_t file, const char* name,
                          const std::vector<double>& values,
                          const std::string& units) {
  const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
  Space space(check_id(H5Screate_simple(1, dims, nullptr), "create dataspace"));
  Dataset dataset(check_id(
      H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space.id(), H5P_DEFAULT,
                 H5P_DEFAULT, H5P_DEFAULT),
      "create dataset " + std::string(name)));
  check_status(H5Dwrite(dataset.id(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, values.data()),
               "write dataset " + std::string(name));
  write_string_attribute(dataset.id(), "units", units);
}

}  // namespace

InputData load_input(const std::filesystem::path& path) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  File file(check_id(
      H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
      "open " + path.string()));

  InputData data;
  data.q1_edges = read_vector(file.id(), "alpha0_edges");
  data.q2_edges = read_vector(file.id(), "logE_edges");
  if (data.q1_edges.size() < 2 || data.q2_edges.size() < 2) {
    throw std::runtime_error("HDF5: coordinate edge arrays are too short");
  }
  data.n1 = data.q1_edges.size() - 1;
  data.n2 = data.q2_edges.size() - 1;
  if (data.n1 > std::numeric_limits<std::size_t>::max() / data.n2) {
    throw std::runtime_error("HDF5: n1*n2 overflows size_t");
  }

  data.q1_centers = read_vector(file.id(), "alpha0_centers");
  data.q2_centers = read_vector(file.id(), "logE_centers");
  if (data.q1_centers.size() != data.n1 ||
      data.q2_centers.size() != data.n2) {
    throw std::runtime_error("HDF5: coordinate center arrays have wrong shape");
  }

  data.G = read_matrix(file.id(), "G", data.n1, data.n2);
  data.D11 = read_matrix(file.id(), "D11", data.n1, data.n2);
  data.D12 = read_matrix(file.id(), "D12", data.n1, data.n2);
  data.D22 = read_matrix(file.id(), "D22", data.n1, data.n2);
  return data;
}

HDF5Output::HDF5Output(const std::filesystem::path& path, std::size_t n1,
                       std::size_t n2, HDF5OutputNames names,
                       HDF5OutputUnits units)
    : n1_(n1), n2_(n2), names_(std::move(names)), units_(std::move(units)) {
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  const hid_t file = check_id(
      H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                H5P_DEFAULT),
      "create " + path.string());
  file_id_ = file;
  try {
    group_id_ = check_id(
        H5Gcreate2(file_id_, names_.fields.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                   H5P_DEFAULT),
        "create group /" + names_.fields);
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

void HDF5Output::write_coordinate1(const std::vector<double>& values) {
  if (values.size() != n1_) {
    throw std::runtime_error("HDF5: /" + names_.coordinate1 +
                             " has the wrong length");
  }
  write_vector_dataset(file_id_, names_.coordinate1.c_str(), values,
                       units_.coordinate1);
}

void HDF5Output::write_coordinate2(const std::vector<double>& values) {
  if (values.size() != n2_) {
    throw std::runtime_error("HDF5: /" + names_.coordinate2 +
                             " has the wrong length");
  }
  write_vector_dataset(file_id_, names_.coordinate2.c_str(), values,
                       units_.coordinate2);
}

void HDF5Output::write_times(const std::vector<double>& values) {
  write_vector_dataset(file_id_, names_.time.c_str(), values,
                       units_.time);
}

void HDF5Output::write_field(std::size_t frame, const Eigen::MatrixXd& field) {
  if (field.rows() != static_cast<Eigen::Index>(n1_) ||
      field.cols() != static_cast<Eigen::Index>(n2_)) {
    throw std::runtime_error("HDF5: /" + names_.fields + "/" +
                             std::to_string(frame) +
                             " has the wrong dimensions");
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
  const std::string name = names_.fields + "/" + std::to_string(frame);
  Dataset dataset(check_id(
      H5Dcreate2(file_id_, name.c_str(), H5T_NATIVE_DOUBLE, space.id(),
                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
      "create dataset /" + name));
  check_status(H5Dwrite(dataset.id(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, row_major.data()),
               "write dataset /" + name);
  write_string_attribute(dataset.id(), "units", units_.fields);
  check_status(H5Fflush(file_id_, H5F_SCOPE_GLOBAL), "flush output file");
}

}  // namespace sayram2d::app
