#pragma once

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace sayram2d {

// Non-owning, strided view of a logical n1 x n2 array.
//
// Strides are measured in elements, not bytes. The logical element (i, j) is
// at data()[i * stride1() + j * stride2()].
template <typename T>
class ArrayView2D {
 public:
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using size_type = std::size_t;
  using stride_type = std::ptrdiff_t;

  ArrayView2D() noexcept = default;

  ArrayView2D(T* data, size_type n1, size_type n2, stride_type s1,
              stride_type s2) noexcept
      : data_(data), n1_(n1), n2_(n2), s1_(s1), s2_(s2) {}

  // Row-major contiguous view.
  ArrayView2D(T* data, size_type n1, size_type n2) noexcept
      : ArrayView2D(data, n1, n2, static_cast<stride_type>(n2), 1) {}

  template <typename U>
  ArrayView2D(const ArrayView2D<U>& other,
              std::enable_if_t<std::is_convertible_v<U*, T*>, int> = 0) noexcept
      : data_(other.data()),
        n1_(other.extent1()),
        n2_(other.extent2()),
        s1_(other.stride1()),
        s2_(other.stride2()) {}

  T* data() const noexcept { return data_; }
  size_type extent1() const noexcept { return n1_; }
  size_type extent2() const noexcept { return n2_; }
  stride_type stride1() const noexcept { return s1_; }
  stride_type stride2() const noexcept { return s2_; }
  size_type size() const noexcept { return n1_ * n2_; }
  bool empty() const noexcept { return n1_ == 0 || n2_ == 0; }

  // True for either contiguous row-major or contiguous column-major storage.
  bool is_contiguous() const noexcept {
    if (empty()) return true;
    return (s1_ == static_cast<stride_type>(n2_) && s2_ == 1) ||
           (s1_ == 1 && s2_ == static_cast<stride_type>(n1_));
  }

  T& operator()(size_type i, size_type j) const noexcept {
    assert(i < n1_ && j < n2_);
    return data_[i * s1_ + j * s2_];
  }

 private:
  T* data_ = nullptr;
  size_type n1_ = 0;
  size_type n2_ = 0;
  stride_type s1_ = 0;
  stride_type s2_ = 0;
};

template <typename T>
ArrayView2D<T> make_view2d(T* data, std::size_t n1, std::size_t n2) {
  return ArrayView2D<T>(data, n1, n2);
}

template <typename T>
ArrayView2D<T> make_view2d(T* data, std::size_t n1, std::size_t n2,
                           std::ptrdiff_t s1, std::ptrdiff_t s2) {
  return ArrayView2D<T>(data, n1, n2, s1, s2);
}

template <typename T>
ArrayView2D<T> make_view2d(std::vector<T>& values, std::size_t n1,
                           std::size_t n2) {
  if (values.size() != n1 * n2) {
    throw std::invalid_argument("make_view2d: vector size does not match extents");
  }
  return ArrayView2D<T>(values.data(), n1, n2);
}

template <typename T>
ArrayView2D<const T> make_view2d(const std::vector<T>& values, std::size_t n1,
                                 std::size_t n2) {
  if (values.size() != n1 * n2) {
    throw std::invalid_argument("make_view2d: vector size does not match extents");
  }
  return ArrayView2D<const T>(values.data(), n1, n2);
}

}  // namespace sayram2d
