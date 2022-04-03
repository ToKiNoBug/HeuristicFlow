/*
 Copyright © 2022  TokiNoBug
This file is part of Heuristic.

    Heuristic is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Heuristic is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Heuristic.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef Heu_MATRIXFIXEDSIZE_H
#define Heu_MATRIXFIXEDSIZE_H

#include <stdint.h>
#include <array>
#include <assert.h>

#include "InternalHeaderCheck.h"

namespace Heu {

template <class Scalar_t, size_t Rows, size_t Cols>
class MatrixFixedSize {
 protected:
  using fast_t = typename std::conditional<sizeof(Scalar_t) <= 3 * sizeof(void *), Scalar_t, const Scalar_t &>::type;

 public:
  MatrixFixedSize(){};
  ~MatrixFixedSize(){};

  using iterator = Scalar_t *;
  using citerator = const Scalar_t *;

  /// Deep copy
  explicit MatrixFixedSize(const MatrixFixedSize &src) {
    for (size_t i = 0; i < size(); i++) {
      array[i] = src.array[i];
    }
  }

  inline iterator begin() noexcept { return array.data(); }

  inline citerator begin() const noexcept { return array.data(); }

  inline iterator end() noexcept { return array.data() + size(); }

  inline citerator end() const noexcept { return array.data() + size(); }

  inline constexpr size_t size() noexcept { return Rows * Cols; }

  inline constexpr size_t rows() noexcept { return Rows; }

  inline constexpr size_t cols() noexcept { return Cols; }

  inline const Scalar_t &operator()(size_t n) const noexcept { return array[n]; }

  inline const Scalar_t &operator[](size_t n) const noexcept { return array[n]; }

  inline const Scalar_t &operator()(size_t r, size_t c) const noexcept { return array[Rows * c + r]; }

  inline Scalar_t &operator()(size_t n) noexcept { return array[n]; }

  inline Scalar_t &operator[](size_t n) noexcept { return array[n]; }

  inline Scalar_t &operator()(size_t r, size_t c) noexcept { return array[Rows * c + r]; }

  inline Scalar_t *data() noexcept { return array.data(); }

  inline const Scalar_t *data() const noexcept { return array.data(); }

  inline static void resize(size_t _r, size_t _c) {
    assert(_r == Rows);
    assert(_c == Cols);
  }

  /// Deep copy
  const MatrixFixedSize &operator=(const MatrixFixedSize &src) {
    for (size_t i = 0; i < size(); i++) {
      array[i] = src.array[i];
    }
    return *this;
  }

  inline void fill(fast_t src) noexcept {
    for (auto &i : *this) {
      i = src;
    }
  }

  static const bool isFixedSize = true;

 protected:
  std::array<Scalar_t, Rows * Cols> array;
};
}  // namespace Heu

#endif  // MATRIXFIXEDSIZE_H
