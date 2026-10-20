/*
 * File:        Mesh.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef MESH_H_
#define MESH_H_

#include <array>
#include <cstddef>
#include <vector>

#include "Edge.h"
#include "Grid2D.h"

struct Ind{
  std::size_t i;
  std::size_t j;
}; 

struct NbrInd {  // neighor index. avoid underflow of i-1 in build_connectivity_
  std::ptrdiff_t i;
  std::ptrdiff_t j;
  bool valid;
};

class Mesh {
  public:

    Mesh(const Grid2D& grid, double dt);

    // Geometry (cell centers)
    const std::vector<double>& x() const { return x_; }
    const std::vector<double>& y() const { return y_; }

    double x(std::size_t i) const { return x_[i]; }
    double y(std::size_t j) const { return y_[j]; }

    std::size_t nx() const { return nx_; }
    std::size_t ny() const { return ny_; }

    // Canonical edges
    double x_edge(std::size_t i) const { return x_edges_[i]; }  // i in [0, nx]
    double y_edge(std::size_t j) const { return y_edges_[j]; }  // j in [0, ny]

    // Nonuniform-safe widths
    double dx(std::size_t i) const { return dx_cell_[i]; }      // i in [0, nx-1]
    double dy(std::size_t j) const { return dy_cell_[j]; }      // j in [0, ny-1]

    double dt() const { return dt_; }

    // double cell_area_dt(std::size_t i, std::size_t j) const {
    //   return dx(i) * dy(j) / dt();
    // }

    double cell_area_dt(const Ind& ind) const {
      return dx(ind.i) * dy(ind.j) / dt();
    }


    std::size_t flatten_cell_index(const Ind& ind) const {
      return ind.j * nx() + ind.i;
    }

    // Neighbor conventions
    // define the neighbor # of 4 adjacent cells
    // im -- (i-1, j, k); jp -- (i, j+1, k)
    // ip -- (i+1, j, k); jm -- (i, j-1, k)
    std::size_t nnbrs() const { return 4; } // each cell has 4 nbrs

    int inbr_im() const { return 0; }
    int inbr_jp() const { return 1; }
    int inbr_ip() const { return 2; }
    int inbr_jm() const { return 3; }

    int rinbr(int inbr) const {
      return rinbr_[static_cast<std::size_t>(inbr)];
    }

    bool get_nbr_ind(const Ind& ind, int inbr, Ind* nbr_indp) const {
      const NbrInd& n = nbr_inds_[idx_(ind.i, ind.j, inbr)];
      if (!n.valid) return false;

      nbr_indp->i = static_cast<std::size_t>(n.i);
      nbr_indp->j = static_cast<std::size_t>(n.j);
      return true;
    }

    void get_nbr_edge(const Ind& ind, int inbr, Edge* edgep) const {
      *edgep = edges_[idx_(ind.i, ind.j, inbr)];
    }

  private:

    std::size_t nx_{0};
    std::size_t ny_{0};
    double dt_{0.0};
  
    // Canonical edges
    std::vector<double> x_edges_;  // size nx+1
    std::vector<double> y_edges_;  // size ny+1

    // Derived geometry
    std::vector<double> dx_cell_;  // size nx
    std::vector<double> dy_cell_;  // size ny
    std::vector<double> x_;        // size nx
    std::vector<double> y_;        // size ny
  
    // Reverse neighbor mapping: im<->ip, jp<->jm
    std::array<int, 4> rinbr_{{2, 3, 0, 1}};

    // Connectivity (flattened): size nx*ny*4
    std::vector<NbrInd> nbr_inds_;
    std::vector<Edge> edges_;

    void build_geometry_from_edges_();
    void build_connectivity_();

    std::size_t cell_id_(std::size_t i, std::size_t j) const {
      return j * nx_ + i;
    }

    std::size_t idx_(std::size_t i, std::size_t j, int inbr) const {
      return 4 * cell_id_(i, j) + static_cast<std::size_t>(inbr);
    }
};

#endif /* MESH_H */

