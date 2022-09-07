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
#include <assert.h>
#include <array>

#include "InternalHeaderCheck.h"

#include "MatrixBase.hpp"

namespace heu {

/**
 * \ingroup HEU_SIMPLEMATRIX
 * \brief Matrix with fixed size
 *
 * \tparam Scalar_t Type of element
 * \tparam Rows Row number
 * \tparam Cols Coloumn number
 */
template <class Scalar, size_t Rows, size_t Cols>
class MatrixFixedSize : public MatrixBase<MatrixFixedSize<Scalar, Rows, Cols>> {
 public:
  using Scalar_t = Scalar;

 protected:
  using fast_t = typename std::conditional<sizeof(Scalar_t) <= 3 * sizeof(void *), Scalar_t,
                                           const Scalar_t &>::type;

 public:
  /**
   * \brief Construct a new Matrix Fixed Size object
   *
   */
  MatrixFixedSize() = default;

  /**
   * \brief Destroy the Matrix Fixed Size object
   *
   */
  ~MatrixFixedSize() = default;

  /// Non-constant iterator (actually a pointer)
  using iterator = Scalar_t *;

  /// Constant iterator (const-ptr actually)
  using citerator = const Scalar_t *;

  /**
   * \brief Deep copy function
   *
   * \param src Source of copying
   */
  explicit MatrixFixedSize(const MatrixFixedSize &src) {
    if constexpr (!std::is_trivially_copy_assignable_v<Scalar_t>) {
      for (size_t i = 0; i < size(); i++) {
        array[i] = src.array[i];
      }
    } else {
      memcpy(array.data(), src.array.data(), sizeof(Scalar_t) * size());
    }
  }

  inline iterator begin() noexcept { return array.data(); }

  inline citerator begin() const noexcept { return array.data(); }

  inline iterator end() noexcept { return array.data() + size(); }

  inline citerator end() const noexcept { return array.data() + size(); }

  constexpr size_t size() const noexcept { return Rows * Cols; }

  constexpr size_t rows() const noexcept { return Rows; }

  constexpr size_t cols() const noexcept { return Cols; }

  inline const Scalar_t &operator()(size_t n) const noexcept { return array[n]; }

  inline const Scalar_t &operator[](size_t n) const noexcept { return array[n]; }

  inline const Scalar_t &operator()(size_t r, size_t c) const noexcept {
    return array[Rows * c + r];
  }

  inline Scalar_t &operator()(size_t n) noexcept { return array[n]; }

  inline Scalar_t &operator[](size_t n) noexcept { return array[n]; }

  inline Scalar_t &operator()(size_t r, size_t c) noexcept { return array[Rows * c + r]; }

  inline Scalar_t *data() noexcept { return array.data(); }

  inline const Scalar_t *data() const noexcept { return array.data(); }

  inline static void resize(size_t _r, size_t _c) noexcept {
    assert(_r == Rows);
    assert(_c == Cols);
  }

  //  operator= for same type
  MatrixFixedSize &operator=(const MatrixFixedSize &src) noexcept {
    for (size_t i = 0; i < size(); i++) {
      array[i] = src.array[i];
    }
    return *this;
  }

  //  operator= for different type
  template <class DerivedB>
  inline MatrixFixedSize &operator=(const MatrixBase<DerivedB> &src) noexcept {
    return MatrixBase<MatrixFixedSize>::template operator=
        <DerivedB>(static_cast<const DerivedB &>(src));
  }

  inline void fill(fast_t src) noexcept {
    for (auto &i : *this) {
      i = src;
    }
  }
  // static constexpr bool isClass = std::is_class<Scalar_t>::value;
  static constexpr bool isFixedSize = true;
  static constexpr int rowsAtCompileTime = Rows;
  static constexpr int colsAtCompileTime = Cols;

 protected:
  std::array<Scalar_t, Rows * Cols> array;
};
}  // namespace heu

#endif  // MATRIXFIXEDSIZE_H
