#pragma once

#include <cstddef>
#include <stdexcept>

#include "sayram2d/array_view2d.h"
#include "sayram2d/boundary_conditions.h"

namespace sayram2d {

// Cell-centered coefficient views plus a copy of the boundary configuration.
// The underlying coefficient and boundary arrays are owned by the caller and
// must remain valid for as long as this Equation is used. Values may be
// updated in place between solves.
class Equation {
 public:
  Equation(ArrayView2D<const double> G, ArrayView2D<const double> D11,
           ArrayView2D<const double> D12, ArrayView2D<const double> D22,
           const BoundaryConditions& bc)
      : G_(G), D11_(D11), D12_(D12), D22_(D22), bc_(bc) {
    if (G_.empty() || D11_.empty() || D12_.empty() || D22_.empty()) {
      throw std::invalid_argument("Equation: coefficient views must be non-empty");
    }
    if (G_.extent1() != D11_.extent1() || G_.extent1() != D12_.extent1() ||
        G_.extent1() != D22_.extent1() || G_.extent2() != D11_.extent2() ||
        G_.extent2() != D12_.extent2() || G_.extent2() != D22_.extent2()) {
      throw std::invalid_argument(
          "Equation: G, D11, D12, and D22 must have identical extents");
    }
  }

  double G(std::size_t i, std::size_t j) const noexcept { return G_(i, j); }
  double D11(std::size_t i, std::size_t j) const noexcept { return D11_(i, j); }
  double D12(std::size_t i, std::size_t j) const noexcept { return D12_(i, j); }
  double D22(std::size_t i, std::size_t j) const noexcept { return D22_(i, j); }

  ArrayView2D<const double> G_view() const noexcept { return G_; }
  ArrayView2D<const double> D11_view() const noexcept { return D11_; }
  ArrayView2D<const double> D12_view() const noexcept { return D12_; }
  ArrayView2D<const double> D22_view() const noexcept { return D22_; }

  void set_loss(ArrayView2D<const double> inv_tau) {
    if (inv_tau.extent1() != G_.extent1() ||
        inv_tau.extent2() != G_.extent2()) {
      throw std::invalid_argument(
          "Equation::set_loss: inverse loss field has the wrong extents");
    }
    inv_tau_ = inv_tau;
    loss_set_ = true;
    loss_constant_ = false;
  }

  void set_loss(double inv_tau) noexcept {
    inv_tau_ = {};
    loss_set_ = true;
    loss_constant_ = true;
    constant_inv_tau_ = inv_tau;
  }

  bool has_loss() const noexcept { return loss_set_; }
  double inv_tau(std::size_t i, std::size_t j) const noexcept {
    if (!loss_set_) return 0.0;
    return loss_constant_ ? constant_inv_tau_ : inv_tau_(i, j);
  }

  ArrayView2D<const double> inv_tau_view() const noexcept {
    return inv_tau_;
  }
  bool loss_is_constant() const noexcept { return loss_constant_; }
  double constant_inv_tau() const noexcept { return constant_inv_tau_; }

  const BoundaryConditions& bc() const noexcept { return bc_; }

 private:
  ArrayView2D<const double> G_;
  ArrayView2D<const double> D11_;
  ArrayView2D<const double> D12_;
  ArrayView2D<const double> D22_;
  ArrayView2D<const double> inv_tau_{};
  double constant_inv_tau_ = 0.0;
  bool loss_set_ = false;
  bool loss_constant_ = true;
  BoundaryConditions bc_;
};

}  // namespace sayram2d
