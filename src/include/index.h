#pragma once

#include <cstddef>

namespace sayram2d {

// Column-major flat index for a contiguous 2D array with n1 rows:
//   (i, j) -> i + j * n1
// This matches Array2D and ArrayView2D.
constexpr std::size_t flat_index(std::size_t i, std::size_t j,
                                 std::size_t n1) noexcept {
  return i + j * n1;
}

}  // namespace sayram2d
