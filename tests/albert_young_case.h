#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "array2d.h"
#include "boundary_conditions.h"
#include "equation.h"
#include "mesh.h"
#include "physics.h"

namespace albert_young {

constexpr double eps = std::numeric_limits<double>::epsilon();

inline double initial_value(double alpha, double log_e) {
  const double energy = std::exp(log_e);
  const double p = sayram2d::physics::momentum_mev_c(energy);
  return std::exp(-(energy - 0.2) / 0.1) *
             (std::sin(alpha) -
              std::sin(sayram2d::physics::degrees_to_radians(5.0))) /
             (p * p) +
         eps;
}

struct Data {
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

inline Data load_data(const std::string& path) {
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
  constexpr char expected_magic[] = "S2DAYD1";
  if (!stream || std::string(magic, 7) != expected_magic || magic[7] != '\0' ||
      version != 1 || n1 <= 0 || n2 <= 0) {
    throw std::runtime_error("invalid Albert-Young binary header: " + path);
  }

  Data data;
  data.n1 = static_cast<std::size_t>(n1);
  data.n2 = static_cast<std::size_t>(n2);
  if (data.n1 > std::numeric_limits<std::size_t>::max() / data.n2) {
    throw std::runtime_error("Albert-Young dimensions overflow: " + path);
  }
  const std::size_t cells = data.n1 * data.n2;

  auto read_array = [&](std::vector<double>& values, std::size_t count,
                        const char* name) {
    values.resize(count);
    stream.read(reinterpret_cast<char*>(values.data()),
                static_cast<std::streamsize>(count * sizeof(double)));
    if (!stream) {
      throw std::runtime_error(std::string("truncated Albert-Young data: ") +
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
  if (stream.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("trailing data after Albert-Young arrays: " + path);
  }
  return data;
}

// Owns the converted cell-centered arrays and the objects whose views point
// into them. No diffusion-coefficient interpolation or unit conversion occurs
// here; that work is done by scripts/convert_albert_young.py.
class Case {
 public:
  explicit Case(const std::string& path)
      : data_(load_data(path)),
        mesh_(data_.q1_edges, data_.q2_edges),
        bc_values_(mesh_.n1() + 1),
        f_(mesh_.n1(), mesh_.n2()) {
    for (std::size_t i = 0; i <= mesh_.n1(); ++i) {
      bc_values_[i] = initial_value(mesh_.q1_edge(i), mesh_.q2_edge(0));
    }

    bc_.set_dirichlet(sayram2d::Side::Q1Min, 0.0);
    bc_.set_zero_flux(sayram2d::Side::Q1Max);
    bc_.set_dirichlet(sayram2d::Side::Q2Min, bc_values_);
    bc_.set_dirichlet(sayram2d::Side::Q2Max, 0.0);

    auto to_array = [&](const std::vector<double>& values) {
      sayram2d::Array2D<double> array(mesh_.n1(), mesh_.n2());
      for (std::size_t i = 0; i < mesh_.n1(); ++i) {
        for (std::size_t j = 0; j < mesh_.n2(); ++j) {
          array(i, j) = values[i * mesh_.n2() + j];
        }
      }
      return array;
    };
    G_ = to_array(data_.G);
    D11_ = to_array(data_.D11);
    D12_ = to_array(data_.D12);
    D22_ = to_array(data_.D22);

    equation_ = std::make_unique<sayram2d::Equation>(
        G_.view(), D11_.view(), D12_.view(), D22_.view(), bc_);

    for (std::size_t i = 0; i < mesh_.n1(); ++i) {
      for (std::size_t j = 0; j < mesh_.n2(); ++j) {
        f_(i, j) = initial_value(mesh_.q1_center(i), mesh_.q2_center(j));
      }
    }
  }

  const sayram2d::Mesh& mesh() const noexcept { return mesh_; }
  sayram2d::Equation& equation() noexcept { return *equation_; }
  const sayram2d::Equation& equation() const noexcept { return *equation_; }
  const Data& data() const noexcept { return data_; }

  sayram2d::Array2D<double>& f() noexcept { return f_; }
  const sayram2d::Array2D<double>& f() const noexcept { return f_; }

 private:
  Data data_;
  sayram2d::Mesh mesh_;
  std::vector<double> bc_values_;
  sayram2d::BoundaryConditions bc_;
  std::unique_ptr<sayram2d::Equation> equation_;
  sayram2d::Array2D<double> G_;
  sayram2d::Array2D<double> D11_;
  sayram2d::Array2D<double> D12_;
  sayram2d::Array2D<double> D22_;
  sayram2d::Array2D<double> f_;
};

}  // namespace albert_young
