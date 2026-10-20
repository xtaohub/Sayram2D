/*
 * File:        Edge.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        01/13/2025 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef EDGE_H_
#define EDGE_H_

#include "common.h"
#include <array>
#include <cstddef>

// Cartesian face normal direction
enum class Direction { XPOS, XNEG, YPOS, YNEG };

using Edge_Vertices = std::array<Point, 2>;

struct VtxInd { std::size_t i, j; };  // vertex index on (x_edges, y_edges)

class Edge{
  public:
    static constexpr int NT = 2;  // number of points for each edge. 

    Edge() : vs_{}, length_(0.0), n_(0.0, 0.0), dir_(Direction::XPOS) {}

    const Point& v(std::size_t i) const { return vs_[i]; } // vertices
    const VtxInd& vind(std::size_t k) const { return vinds_[k]; }

    std::size_t v_size() const { return vs_.size(); }

    double length() const { return length_; }
    const Vector2& n() const { return n_; } // normal vector, depending on the edge direction
    Direction dir() const { return dir_; }

    void set(const Edge_Vertices& vs, Direction dir, const std::array<VtxInd,2>& vinds) {
      vs_ = vs;
      vinds_ = vinds;
      length_ = (vs_[1] - vs_[0]).norm();
      dir_ = dir;

      switch (dir_) {
        case Direction::XPOS: n_ = Vector2( 1.0,  0.0); break;
        case Direction::XNEG: n_ = Vector2(-1.0,  0.0); break;
        case Direction::YPOS: n_ = Vector2( 0.0,  1.0); break;
        case Direction::YNEG: n_ = Vector2( 0.0, -1.0); break;
        default:              n_ = Vector2( 0.0,  0.0); break;  // defensive
      }
    }

  private:

    Edge_Vertices vs_; 
    std::array<VtxInd,2> vinds_;

    double length_; 
    Vector2 n_; // the unit normal vector
    Direction dir_; 
};

#endif /* EDGE_H */

