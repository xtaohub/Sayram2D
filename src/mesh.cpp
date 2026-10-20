#include "mesh.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace sayram2d {
namespace {

void validate_edges(const std::vector<double>& edges, const char* name) {
  if (edges.size() < 2) {
    throw std::invalid_argument(std::string(name) +
                                " must contain at least two edges");
  }
  for (std::size_t k = 0; k < edges.size(); ++k) {
    if (!std::isfinite(edges[k])) {
      throw std::invalid_argument(std::string(name) +
                                  " contains a non-finite coordinate");
    }
    if (k > 0 && !(edges[k] > edges[k - 1])) {
      throw std::invalid_argument(std::string(name) +
                                  " must be strictly increasing");
    }
  }
}

void fill_centers_and_widths(const std::vector<double>& edges,
                             std::vector<double>* centers,
                             std::vector<double>* widths) {
  const std::size_t n = edges.size() - 1;
  centers->resize(n);
  widths->resize(n);
  for (std::size_t k = 0; k < n; ++k) {
    (*centers)[k] = 0.5 * (edges[k] + edges[k + 1]);
    (*widths)[k] = edges[k + 1] - edges[k];
  }
}

}  // namespace

Mesh::Mesh(std::vector<double> q1_edges, std::vector<double> q2_edges)
    : q1_edges_(std::move(q1_edges)),
      q2_edges_(std::move(q2_edges)) {
  validate_edges(q1_edges_, "q1_edges");
  validate_edges(q2_edges_, "q2_edges");

  n1_ = q1_edges_.size() - 1;
  n2_ = q2_edges_.size() - 1;

  fill_centers_and_widths(q1_edges_, &q1_centers_, &dq1_);
  fill_centers_and_widths(q2_edges_, &q2_centers_, &dq2_);
}

bool Mesh::neighbor(Ind cell, Side side, Ind* result) const noexcept {
  assert(cell.i < n1_ && cell.j < n2_);
  switch (side) {
    case Side::Q1Min:
      if (cell.i == 0) return false;
      *result = {cell.i - 1, cell.j};
      return true;
    case Side::Q1Max:
      if (cell.i + 1 >= n1_) return false;
      *result = {cell.i + 1, cell.j};
      return true;
    case Side::Q2Min:
      if (cell.j == 0) return false;
      *result = {cell.i, cell.j - 1};
      return true;
    case Side::Q2Max:
      if (cell.j + 1 >= n2_) return false;
      *result = {cell.i, cell.j + 1};
      return true;
  }
  return false;
}

Face Mesh::face(Ind cell, Side side) const {
  assert(cell.i < n1_ && cell.j < n2_);

  Face result;
  result.owner = cell;
  result.side = side;

  const std::size_t i = cell.i;
  const std::size_t j = cell.j;
  const double q1l = q1_edges_[i];
  const double q1r = q1_edges_[i + 1];
  const double q2l = q2_edges_[j];
  const double q2r = q2_edges_[j + 1];

  switch (side) {
    case Side::Q1Min:
      result.a = {q1l, q2r};
      result.b = {q1l, q2l};
      result.va = {i, j + 1};
      result.vb = {i, j};
      result.length = dq2_[j];
      result.normal = {-1.0, 0.0};
      if (i > 0) {
        result.boundary = false;
        result.neighbor = {i - 1, j};
      }
      break;
    case Side::Q1Max:
      result.a = {q1r, q2l};
      result.b = {q1r, q2r};
      result.va = {i + 1, j};
      result.vb = {i + 1, j + 1};
      result.length = dq2_[j];
      result.normal = {1.0, 0.0};
      if (i + 1 < n1_) {
        result.boundary = false;
        result.neighbor = {i + 1, j};
      }
      break;
    case Side::Q2Min:
      result.a = {q1l, q2l};
      result.b = {q1r, q2l};
      result.va = {i, j};
      result.vb = {i + 1, j};
      result.length = dq1_[i];
      result.normal = {0.0, -1.0};
      if (j > 0) {
        result.boundary = false;
        result.neighbor = {i, j - 1};
      }
      break;
    case Side::Q2Max:
      result.a = {q1r, q2r};
      result.b = {q1l, q2r};
      result.va = {i + 1, j + 1};
      result.vb = {i, j + 1};
      result.length = dq1_[i];
      result.normal = {0.0, 1.0};
      if (j + 1 < n2_) {
        result.boundary = false;
        result.neighbor = {i, j + 1};
      }
      break;
  }
  return result;
}

}  // namespace sayram2d
