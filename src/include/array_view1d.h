#pragma once

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace sayram2d {

// Non-owning, strided view of a logical 1D array.
//
// Strides are measured in elements, not bytes. The logical element (i) is at
// data()[i * stride()].
template <typename T>
class ArrayView1D {
 public:
  using element_type = T;
  using value_type = std::remove_cv_t<T>;
  using size_type = std::size_t;
  using stride_type = std::ptrdiff_t;

  ArrayView1D() noexcept = default;

  ArrayView1D(T* data, size_type extent, stride_type stride) noexcept
      : data_(data), extent_(extent), stride_(stride) {}

  // Contiguous view.
  ArrayView1D(T* data, size_type extent) noexcept
      : ArrayView1D(data, extent, 1) {}

  template <typename U>
  ArrayView1D(const ArrayView1D<U>& other,
              std::enable_if_t<std::is_convertible_v<U*, T*>, int> = 0) noexcept
      : data_(other.data()),
        extent_(other.extent()),
        stride_(other.stride()) {}

  T* data() const noexcept { return data_; }
  size_type extent() const noexcept { return extent_; }
  size_type size() const noexcept { return extent_; }
  stride_type stride() const noexcept { return stride_; }
  bool empty() const noexcept { return extent_ == 0; }

  bool is_contiguous() const noexcept {
    return empty() || stride_ == 1;
  }

  T& operator()(size_type i) const noexcept {
    assert(i < extent_);
    return data_[i * stride_];
  }

 private:
  T* data_ = nullptr;
  size_type extent_ = 0;
  stride_type stride_ = 0;
};

template <typename T>
ArrayView1D<T> make_view1d(T* data, std::size_t extent,
                           std::ptrdiff_t stride) {
  return ArrayView1D<T>(data, extent, stride);
}

template <typename T>
ArrayView1D<T> make_view1d(T* data, std::size_t extent) {
  return ArrayView1D<T>(data, extent);
}

}  // namespace sayram2d
