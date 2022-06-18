/*
 Copyright © 2022  TokiNoBug
This file is part of HeuristicFlow.

    HeuristicFlow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HeuristicFlow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HeuristicFlow.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef Heu_MATRIXDYNAMICSIZE_H
#define Heu_MATRIXDYNAMICSIZE_H

#include <stdint.h>
#include <memory>
#include <type_traits>

#include "InternalHeaderCheck.h"

namespace heu {
/// col-major matrix with basic types only

/**
 * \ingroup HEU_SIMPLEMATRIX
 * \brief Matrix with dynamic size
 * \class MatrixDynamicSize
 *
 * This template class maintains a matrix with dynamic size. It supports all kinds of element but it doesn't implement
 * any math computations.
 *
 * It has both STL-like APIs and Eigen-like APIs.
 *
 * \tparam Scalar_t Type of element
 * \tparam allocator_t Type of allocator
 */
template <class Scalar_t, class allocator_t = std::allocator<Scalar_t>>
class MatrixDynamicSize {
 protected:
  using fast_t = typename std::conditional<sizeof(Scalar_t) <= 3 * sizeof(void *), Scalar_t, const Scalar_t &>::type;

 public:
  /**
   * \brief Construct a new Matrix Dynamic Size object. Set the size to 0 and dataPtr to nullptr
   *
   */
  MatrixDynamicSize() {
    rowNum = 0;
    colNum = 0;
    _capacity = 0;
    dataPtr = nullptr;
  };

  /**
   * \brief Construct a new Matrix Dynamic Size object with given rows and cols.
   *
   * \param r
   * \param c
   */
  MatrixDynamicSize(size_t r, size_t c) {
    rowNum = r;
    colNum = c;
    _capacity = r * c;
    dataPtr = alloc().allocate(_capacity);
    if (isClass)
      for (size_t i = 0; i < _capacity; i++) {
        alloc().construct(dataPtr + i);
      }
  }

  /**
   * \brief Deep copy constructor
   *
   * \param src Source of copying
   */
  explicit MatrixDynamicSize(const MatrixDynamicSize &src) {
    resize(src.rowNum, src.colNum);
    for (size_t i = 0; i < size(); i++) {
      dataPtr[i] = src.dataPtr[i];
    }
  }

  /// Deep copy

  /**
   * \brief Deep copy operator
   *
   * \param src Source of copying
   * \return const MatrixDynamicSize& A const ref to *this
   */
  const MatrixDynamicSize &operator=(const MatrixDynamicSize &src) {
    resize(src.rowNum, src.colNum);
    for (size_t i = 0; i < size(); i++) {
      dataPtr[i] = src.dataPtr[i];
    }
    return *this;
  }

  /**
   * \brief Destroy the Matrix Dynamic Size object and all its elements
   *
   */
  ~MatrixDynamicSize() {
    if (dataPtr != nullptr) {
      alloc().deallocate(dataPtr, _capacity);
      if (isClass)
        for (size_t i = 0; i < _capacity; i++) {
          alloc().destroy(dataPtr + i);
        }
    }
  };

  /**
   * \brief Type of non-const iterator
   *
   */
  using iterator = Scalar_t *;

  /**
   * \brief Type of const iterator
   *
   */
  using citerator = const Scalar_t *;

  /**
   * \brief Same as that in STL
   *
   */
  inline iterator begin() noexcept { return dataPtr; }

  /**
   * \brief Same as that in STL
   *
   */
  inline iterator end() noexcept { return dataPtr + size(); }

  /**
   * \brief Same as that in STL
   *
   */
  inline citerator begin() const noexcept { return dataPtr; }

  /**
   * \brief Same as that in STL
   *
   */
  inline citerator end() const noexcept { return dataPtr + size(); }

  /**
   * \brief Number of elements. Same as that in STL
   *
   */
  inline size_t size() const noexcept { return rowNum * colNum; }

  /**
   * \brief To reserve some space. Same as that in std::vector
   *
   * This function will do nothing if s is less than current capacity
   *
   * \param s Size to be reserved
   */
  void reserve(size_t s) {
    if (_capacity <= s) return;

    if (isClass)
      for (size_t i = 0; i < _capacity; i++) {
        alloc().destroy(dataPtr + i);
      }
    alloc().deallocate(dataPtr, _capacity);

    dataPtr = alloc().allocate(s);
    _capacity = s;
    if (isClass) {
      for (size_t i = 0; i < _capacity; i++) {
        alloc().construct(dataPtr + i);
      }
    }
  }

  /**
   * \brief Change the shape of matrix
   *
   * \param r New row numbers
   * \param c New coloumn numbers
   */
  void resize(size_t r, size_t c) {
    if (r * c != size()) {
      if (r * c > _capacity) {
        if (dataPtr != nullptr) {
          if (isClass)
            for (size_t i = 0; i < _capacity; i++) {
              alloc().destroy(dataPtr + i);
            }
          alloc().deallocate(dataPtr, _capacity);
        }

        dataPtr = alloc().allocate(r * c);
        _capacity = r * c;

        if (isClass)
          for (size_t i = 0; i < _capacity; i++) {
            alloc().construct(dataPtr + i);
          }
      }
    }

    rowNum = r;
    colNum = c;
  }

  /**
   * \brief Same API to eigen's matrices
   */
  inline size_t rows() const noexcept { return rowNum; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline size_t cols() const noexcept { return colNum; }

  /**
   * \brief Same API to std::vector
   *
   */
  inline size_t capacity() const noexcept { return _capacity; }

  /**
   * \brief Same API to std::vector
   */
  inline const Scalar_t &operator[](size_t n) const noexcept { return dataPtr[n]; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline const Scalar_t &operator()(size_t n) const noexcept { return dataPtr[n]; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline const Scalar_t &operator()(size_t r, size_t c) const noexcept { return dataPtr[rowNum * c + r]; }

  /**
   * \brief Same API to std::vector
   */
  inline Scalar_t &operator[](size_t n) noexcept { return dataPtr[n]; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline Scalar_t &operator()(size_t n) noexcept { return dataPtr[n]; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline Scalar_t &operator()(size_t r, size_t c) noexcept { return dataPtr[rowNum * c + r]; }

  /**
   * \brief Same API to std::vector
   */
  inline Scalar_t *data() noexcept { return dataPtr; }

  /**
   * \brief Same API to std::vector
   */
  inline const Scalar_t *data() const noexcept { return dataPtr; }

  /**
   * \brief Fill the matrix with a element
   *
   * \param src
   */
  inline void fill(fast_t src) noexcept {
    for (auto &i : *this) {
      i = src;
    }
  }

  static const bool isFixedSize = false;

 protected:
  Scalar_t *dataPtr;
  size_t rowNum;
  size_t colNum;
  size_t _capacity;

 private:
  static constexpr bool isClass = std::is_class<Scalar_t>::value;

  inline static allocator_t &alloc() {
    static allocator_t alloctor;
    return alloctor;
  }
};

}  // namespace heu

#endif  // MATRIXDYNAMICSIZE_H
