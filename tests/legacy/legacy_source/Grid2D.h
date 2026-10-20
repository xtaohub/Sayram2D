/*
 * File:        Grid2D.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        12/30/2025 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#pragma once

#include <vector>
#include <stdexcept>
#include <string>

class Grid2D {
public:
  // Edge coordinates (canonical representation)
  // x_edges: size nx+1
  // y_edges: size ny+1
  std::vector<double> x_edges;
  std::vector<double> y_edges;

  Grid2D() = default;

  Grid2D(std::vector<double> x_edges_in,
         std::vector<double> y_edges_in)
    : x_edges(std::move(x_edges_in)),
      y_edges(std::move(y_edges_in)) {
    validate();
  }

  // number of cells
  std::size_t nx() const { return x_edges.size() - 1; }
  std::size_t ny() const { return y_edges.size() - 1; }

  // domain bounds (Cartesian boundaries)
  double x_min() const { return x_edges.front(); }
  double x_max() const { return x_edges.back(); }
  double y_min() const { return y_edges.front(); }
  double y_max() const { return y_edges.back(); }

  // basic sanity checks
  void validate() const {
    if (x_edges.size() < 2) {
      throw std::runtime_error("Grid2D: x_edges must have size >= 2.");
    }
    if (y_edges.size() < 2) {
      throw std::runtime_error("Grid2D: y_edges must have size >= 2.");
    }

    for (std::size_t i = 0; i + 1 < x_edges.size(); ++i) {
      if (!(x_edges[i + 1] > x_edges[i])) {
        throw std::runtime_error(
          "Grid2D: x_edges must be strictly increasing at i=" +
          std::to_string(i));
      }
    }

    for (std::size_t j = 0; j + 1 < y_edges.size(); ++j) {
      if (!(y_edges[j + 1] > y_edges[j])) {
        throw std::runtime_error(
          "Grid2D: y_edges must be strictly increasing at j=" +
          std::to_string(j));
      }
    }
  }
};


