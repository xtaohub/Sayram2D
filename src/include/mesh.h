#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "index.h"

namespace sayram2d {

enum class Side : std::uint8_t { Q1Min = 0, Q1Max = 1, Q2Min = 2, Q2Max = 3 };

constexpr Side opposite(Side side) noexcept {
  switch (side) {
    case Side::Q1Min: return Side::Q1Max;
    case Side::Q1Max: return Side::Q1Min;
    case Side::Q2Min: return Side::Q2Max;
    case Side::Q2Max: return Side::Q2Min;
  }
  return Side::Q1Min;
}

constexpr std::array<Side, 4> all_sides() noexcept {
  return {Side::Q1Min, Side::Q1Max, Side::Q2Min, Side::Q2Max};
}

struct Ind {
  std::size_t i = 0;
  std::size_t j = 0;
};

struct VtxInd {
  std::size_t i = 0;
  std::size_t j = 0;
};

struct Vec2 {
  double q1 = 0.0;
  double q2 = 0.0;
};

struct Face {
  Ind owner{};
  Side side = Side::Q1Min;
  bool boundary = true;
  Ind neighbor{};
  VtxInd va{};
  VtxInd vb{};
  Vec2 a{};
  Vec2 b{};
  double length = 0.0;
  Vec2 normal{};
};

class Mesh {
 public:
  Mesh(std::vector<double> q1_edges, std::vector<double> q2_edges);

  std::size_t n1() const noexcept { return n1_; }
  std::size_t n2() const noexcept { return n2_; }
  std::size_t n_cells() const noexcept { return n1_ * n2_; }
  std::size_t n_vertices1() const noexcept { return n1_ + 1; }
  std::size_t n_vertices2() const noexcept { return n2_ + 1; }
  std::size_t cell_index(Ind cell) const noexcept {
    return flat_index(cell.i, cell.j, n2_);
  }

  const std::vector<double>& q1_edges() const noexcept { return q1_edges_; }
  const std::vector<double>& q2_edges() const noexcept { return q2_edges_; }
  const std::vector<double>& q1_centers() const noexcept { return q1_centers_; }
  const std::vector<double>& q2_centers() const noexcept { return q2_centers_; }

  double q1_edge(std::size_t i) const noexcept { return q1_edges_[i]; }
  double q2_edge(std::size_t j) const noexcept { return q2_edges_[j]; }
  double q1_center(std::size_t i) const noexcept { return q1_centers_[i]; }
  double q2_center(std::size_t j) const noexcept { return q2_centers_[j]; }
  double dq1(std::size_t i) const noexcept { return dq1_[i]; }
  double dq2(std::size_t j) const noexcept { return dq2_[j]; }

  Vec2 cell_center(Ind cell) const noexcept {
    return {q1_centers_[cell.i], q2_centers_[cell.j]};
  }
  double cell_measure(Ind cell) const noexcept {
    return dq1_[cell.i] * dq2_[cell.j];
  }
  Vec2 vertex_position(VtxInd vertex) const noexcept {
    return {q1_edges_[vertex.i], q2_edges_[vertex.j]};
  }

  // Returns false at a boundary. `neighbor` is unchanged in that case.
  bool neighbor(Ind cell, Side side, Ind* neighbor) const noexcept;

  // Owner-oriented face geometry. Precondition: valid cell and side.
  // Endpoints are ordered so that clockwise rotation of (b - a) is the
  // outward normal. For the neighbor-oriented face, use
  // face(*neighbor, opposite(side)).
  Face face(Ind cell, Side side) const;

  bool is_boundary_face(Ind cell, Side side) const noexcept {
    Ind ignored{};
    return !neighbor(cell, side, &ignored);
  }
  bool is_boundary_vertex(VtxInd vertex) const noexcept {
    return vertex.i == 0 || vertex.i == n1_ || vertex.j == 0 ||
           vertex.j == n2_;
  }

 private:
  std::size_t n1_ = 0;
  std::size_t n2_ = 0;
  std::vector<double> q1_edges_;
  std::vector<double> q2_edges_;
  std::vector<double> q1_centers_;
  std::vector<double> q2_centers_;
  std::vector<double> dq1_;
  std::vector<double> dq2_;
};

}  // namespace sayram2d
