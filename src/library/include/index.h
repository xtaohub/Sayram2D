#pragma once

#include <cstddef>

namespace sayram2d {

// Flat index for contiguous internal data whose first dimension has extent n1
// and whose second dimension is faster-varying in memory. This is the layout
// used by Mesh cell indexing and the Solver vertex/workspace arrays.
constexpr std::size_t flat_index(std::size_t i, std::size_t j,
                                 std::size_t n1) noexcept {
  return j * n1 + i;
}

}  // namespace sayram2d
