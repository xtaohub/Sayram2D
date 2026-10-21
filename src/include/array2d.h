#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

#include "array_view2d.h"

namespace sayram2d {

// Owning, contiguous, column-major 2D array.
//
// Array2D owns its storage and always uses column-major indexing:
//   value(i, j) == storage[i + j * n1]
// Use view() to pass it to Equation or Solver.
template <typename T>
class Array2D {
 public:
  using value_type = T;
  using size_type = std::size_t;

  Array2D() = default;

  Array2D(size_type n1, size_type n2, T value = T{})
      : n1_(n1), n2_(n2), values_(n1 * n2, value) {}

  size_type extent1() const noexcept { return n1_; }
  size_type extent2() const noexcept { return n2_; }
  size_type size() const noexcept { return values_.size(); }
  bool empty() const noexcept { return values_.empty(); }

  T& operator()(size_type i, size_type j) noexcept {
    assert(i < n1_ && j < n2_);
    return values_[i + j * n1_];
  }

  const T& operator()(size_type i, size_type j) const noexcept {
    assert(i < n1_ && j < n2_);
    return values_[i + j * n1_];
  }

  void fill(const T& value) { std::fill(values_.begin(), values_.end(), value); }

  auto begin() noexcept { return values_.begin(); }
  auto end() noexcept { return values_.end(); }
  auto begin() const noexcept { return values_.begin(); }
  auto end() const noexcept { return values_.end(); }

  ArrayView2D<T> view() & noexcept {
    return ArrayView2D<T>(values_.data(), n1_, n2_);
  }

  ArrayView2D<const T> view() const & noexcept {
    return ArrayView2D<const T>(values_.data(), n1_, n2_);
  }

  ArrayView2D<T> view() && = delete;
  ArrayView2D<const T> view() const && = delete;

 private:
  size_type n1_ = 0;
  size_type n2_ = 0;
  std::vector<T> values_;
};

}  // namespace sayram2d
