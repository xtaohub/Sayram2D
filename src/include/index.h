#pragma once

#include <cstddef>

namespace sayram2d {

// Row-major flat index for a contiguous 2D array with n2 columns:
//   (i, j) -> i * n2 + j
// This matches Array2D and the HDF5 field layout.
constexpr std::size_t flat_index(std::size_t i, std::size_t j,
                                 std::size_t n2) noexcept {
  return i * n2 + j;
}

}  // namespace sayram2d
