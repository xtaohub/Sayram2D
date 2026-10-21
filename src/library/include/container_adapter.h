#pragma once

#include <cstddef>
#include <stdexcept>
#include <type_traits>

#include "array_view2d.h"

namespace sayram2d {

// Adapter for xtensor-like 2-D containers (for example xt::xtensor or
// xt::xarray). The core does not include xtensor; this header is opt-in.
// Requires: data(), dimension(), shape(k), and strides() with operator[].
template <typename Container>
auto make_view2d(Container& container)
    -> ArrayView2D<std::remove_reference_t<decltype(*container.data())>> {
  if (container.dimension() != 2) {
    throw std::invalid_argument("make_view2d: container must be two-dimensional");
  }
  using T = std::remove_reference_t<decltype(*container.data())>;
  const auto strides = container.strides();
  return ArrayView2D<T>(container.data(), container.shape(0),
                        container.shape(1), strides[0], strides[1]);
}

}  // namespace sayram2d
