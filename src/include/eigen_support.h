#pragma once

#include <cstddef>
#include <type_traits>

#include <Eigen/Dense>

#include "array_view2d.h"

namespace sayram2d {

// Zero-copy views of ordinary column-major Eigen::MatrixXd objects.
// n1 = rows, n2 = cols, stride1 = 1, stride2 = n1.
inline ArrayView2D<double> make_view2d(Eigen::MatrixXd& matrix) noexcept {
  return ArrayView2D<double>(matrix.data(),
                             static_cast<std::size_t>(matrix.rows()),
                             static_cast<std::size_t>(matrix.cols()), 1,
                             static_cast<std::ptrdiff_t>(matrix.rows()));
}

inline ArrayView2D<const double> make_view2d(
    const Eigen::MatrixXd& matrix) noexcept {
  return ArrayView2D<const double>(matrix.data(),
                                   static_cast<std::size_t>(matrix.rows()),
                                   static_cast<std::size_t>(matrix.cols()), 1,
                                   static_cast<std::ptrdiff_t>(matrix.rows()));
}

// Reject binding a view to a temporary matrix.
ArrayView2D<const double> make_view2d(Eigen::MatrixXd&&) = delete;

template <typename T>
inline constexpr bool is_eigen_matrix_lvalue_v =
    std::is_lvalue_reference_v<T> &&
    std::is_same_v<std::remove_cv_t<std::remove_reference_t<T>>,
                   Eigen::MatrixXd>;

}  // namespace sayram2d
