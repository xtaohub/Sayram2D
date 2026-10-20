/*
 * File:        Mesh.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        07/19/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Mesh.h"
// #include "Edge.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

Mesh::Mesh(const Grid2D& grid, double dt)
  : nx_(grid.nx()),
    ny_(grid.ny()),
    dt_(dt),
    x_edges_(grid.x_edges),
    y_edges_(grid.y_edges),
    dx_cell_(nx_),
    dy_cell_(ny_),
    x_(nx_),
    y_(ny_) {

  build_geometry_from_edges_();

  nbr_inds_.resize(nx_ * ny_ * nnbrs());
  edges_.resize(nx_ * ny_ * nnbrs());

  build_connectivity_();
}

void Mesh::build_geometry_from_edges_() {
  if (x_edges_.size() != nx_ + 1) {
    throw std::runtime_error("Mesh: x_edges size mismatch.");
  }
  if (y_edges_.size() != ny_ + 1) {
    throw std::runtime_error("Mesh: y_edges size mismatch.");
  }

  for (std::size_t i = 0; i < nx_; ++i) {
    const double xl = x_edges_[i];
    const double xr = x_edges_[i + 1];
    const double dx = xr - xl;
    if (!(dx > 0.0)) {
      throw std::runtime_error("Mesh: non-positive dx at i=" + std::to_string(i));
    }
    dx_cell_[i] = dx;
    x_[i] = 0.5 * (xl + xr);
  }

  for (std::size_t j = 0; j < ny_; ++j) {
    const double yb = y_edges_[j];
    const double yt = y_edges_[j + 1];
    const double dy = yt - yb;
    if (!(dy > 0.0)) {
      throw std::runtime_error("Mesh: non-positive dy at j=" + std::to_string(j));
    }
    dy_cell_[j] = dy;
    y_[j] = 0.5 * (yb + yt);
  }
}

void Mesh::build_connectivity_() {
  // Neighbor ordering:
  // 0: im (west), 1: jp (north), 2: ip (east), 3: jm (south)
  //
  
  auto valid = [](std::ptrdiff_t ii, std::ptrdiff_t jj) -> NbrInd {
    return NbrInd{ii, jj, true};
  };
  auto invalid = []() -> NbrInd {
    return NbrInd{-1, -1, false};
  };

  for (std::size_t i = 0; i < nx_; ++i) {
    for (std::size_t j = 0; j < ny_; ++j) {

      // nbr 0: im (west), XNEG
      {
        const int inbr = 0;

        if (i == 0) {
          nbr_inds_[idx_(i, j, inbr)] = invalid();
        } else {
          nbr_inds_[idx_(i, j, inbr)] =
              valid(static_cast<std::ptrdiff_t>(i) - 1,
                    static_cast<std::ptrdiff_t>(j));
        }

        const Point A = {x_edge(i), y_edge(j + 1)};
        const Point B = {x_edge(i), y_edge(j)};
        edges_[idx_(i, j, inbr)].set({A, B}, Direction::XNEG,
                                    {{{i, j + 1}, {i, j}}});
      }

      // nbr 1: jp (north), YPOS
      {
        const int inbr = 1;

        if (j + 1 >= ny_) {
          nbr_inds_[idx_(i, j, inbr)] = invalid();
        } else {
          nbr_inds_[idx_(i, j, inbr)] =
              valid(static_cast<std::ptrdiff_t>(i),
                    static_cast<std::ptrdiff_t>(j) + 1);
        }

        const Point B = {x_edge(i),     y_edge(j + 1)};
        const Point A = {x_edge(i + 1), y_edge(j + 1)};
        edges_[idx_(i, j, inbr)].set({A, B}, Direction::YPOS,
                                    {{{i + 1, j + 1}, {i, j + 1}}});
      }

      // nbr 2: ip (east), XPOS
      {
        const int inbr = 2;

        if (i + 1 >= nx_) {
          nbr_inds_[idx_(i, j, inbr)] = invalid();
        } else {
          nbr_inds_[idx_(i, j, inbr)] =
              valid(static_cast<std::ptrdiff_t>(i) + 1,
                    static_cast<std::ptrdiff_t>(j));
        }

        const Point B = {x_edge(i + 1), y_edge(j + 1)};
        const Point A = {x_edge(i + 1), y_edge(j)};
        edges_[idx_(i, j, inbr)].set({A, B}, Direction::XPOS,
                                    {{{i + 1, j}, {i + 1, j + 1}}});
      }

      // nbr 3: jm (south), YNEG
      {
        const int inbr = 3;

        if (j == 0) {
          nbr_inds_[idx_(i, j, inbr)] = invalid();
        } else {
          nbr_inds_[idx_(i, j, inbr)] =
              valid(static_cast<std::ptrdiff_t>(i),
                    static_cast<std::ptrdiff_t>(j) - 1);
        }

        const Point B = {x_edge(i + 1), y_edge(j)};
        const Point A = {x_edge(i),     y_edge(j)};
        edges_[idx_(i, j, inbr)].set({A, B}, Direction::YNEG,
                                    {{{i, j}, {i + 1, j}}});
      }
    }
  }
}
