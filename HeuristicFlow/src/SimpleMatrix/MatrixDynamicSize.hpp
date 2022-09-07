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

#include "MatrixBase.hpp"

#include <Eigen/Dense>

namespace heu {
/// col-major matrix with basic types only

/**
 * \ingroup HEU_SIMPLEMATRIX
 * \brief Matrix with dynamic size
 * \class MatrixDynamicSize
 *
 * This template class maintains a matrix with dynamic size. It supports all kinds of element but it
 * doesn't implement any math computations.
 *
 * It has both STL-like APIs and Eigen-like APIs.
 *
 * \tparam Scalar_t Type of element
 * \tparam allocator_t Type of allocator
 */
template <class Scalar, class allocator_t = std::allocator<Scalar>>
class MatrixDynamicSize : public MatrixBase<MatrixDynamicSize<Scalar, allocator_t>> {
 public:
  using Scalar_t = Scalar;

 protected:
  using fast_t = typename std::conditional<sizeof(Scalar_t) <= 3 * sizeof(void *), Scalar_t,
                                           const Scalar_t &>::type;

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
  MatrixDynamicSize(int r, int c) noexcept {
    rowNum = r;
    colNum = c;
    _capacity = r * c;
    dataPtr = alloc().allocate(_capacity);
    if constexpr (!std::is_trivially_default_constructible_v<Scalar_t>)
      for (int i = 0; i < _capacity; i++) {
        // alloc().construct(dataPtr + i);
        dataPtr[i] = Scalar_t();
      }
  }

  /**
   * \brief Deep copy constructor
   *
   * \param src Source of copying
   */
  explicit MatrixDynamicSize(const MatrixDynamicSize &src) noexcept {
    resize(src.rowNum, src.colNum);
    for (int i = 0; i < size(); i++) {
      dataPtr[i] = src.dataPtr[i];
    }
  }

  explicit MatrixDynamicSize(MatrixDynamicSize &&src) noexcept {
    dataPtr = src.dataPtr;
    rowNum = src.rowNum;
    colNum = src.colNum;
    _capacity = src._capacity;

    src.dataPtr = nullptr;
    src.rowNum = 0;
    src.colNum = 0;
    src._capacity = 0;
  }

  /**
   * \brief Destroy the Matrix Dynamic Size object and all its elements
   *
   */
  ~MatrixDynamicSize() {
    if (dataPtr != nullptr) {
      if constexpr (!std::is_trivially_destructible_v<Scalar_t>)
        for (int i = 0; i < _capacity; i++) {
          (dataPtr + i)->~Scalar_t();
          // alloc().destroy(dataPtr + i);
        }

      alloc().deallocate(dataPtr, _capacity);
    }
  };

  /// Deep copy

  /**
   * \brief Deep copy operator
   *
   * \param src Source of copying
   * \return const MatrixDynamicSize& A const ref to *this
   */
  MatrixDynamicSize &operator=(const MatrixDynamicSize &src) noexcept {
    resize(src.rowNum, src.colNum);
    for (int i = 0; i < size(); i++) {
      if constexpr (!std::is_trivially_copy_assignable_v<Scalar_t>) {
        dataPtr[i] = src.dataPtr[i];
      } else {
        memcpy(dataPtr, src, sizeof(Scalar_t) * src.size());
      }
    }
    return *this;
  }

  /**
   * \brief Move assigning operator
   *
   * \param src A rvalue reference
   * \return MatrixDynamicSize&
   */
  MatrixDynamicSize &operator=(MatrixDynamicSize &&src) noexcept {
    if (this->dataPtr != nullptr) {
      if constexpr (!std::is_trivially_destructible_v<Scalar_t>) {
        for (int i = 0; i < _capacity; i++) {
          (dataPtr + i)->~Scalar_t();
          // alloc().destroy(dataPtr + i);
        }
      }
      alloc().deallocate(dataPtr, _capacity);
    }

    this->dataPtr = src.dataPtr;
    src.dataPtr = nullptr;
    this->rowNum = src.rowNum;
    src.rowNum = 0;
    this->colNum = src.colNum;
    src.colNum = 0;
    this->_capacity = src._capacity;
    src._capacity = 0;
  }

  template <class DerivedB>
  inline MatrixDynamicSize &operator=(const MatrixBase<DerivedB> &src) noexcept {
    return MatrixBase<MatrixDynamicSize>::template operator=
        <DerivedB>(static_cast<const DerivedB &>(src));
  }

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
  inline int size() const noexcept { return rowNum * colNum; }

  /**
   * \brief To reserve some space. Same as that in std::vector
   *
   * This function will do nothing if s is less than current capacity
   *
   * \param s Size to be reserved
   */
  void reserve(int s) noexcept {
    if (_capacity <= s) return;

    Scalar_t *newPtr = alloc().allocate(s);

    memcpy(newPtr, dataPtr, sizeof(Scalar_t) * _capacity);

    alloc().deallocate(dataPtr, _capacity);

    _capacity = s;

    dataPtr = newPtr;
  }

  /**
   * \brief Change the shape of matrix
   *
   * \param r New row numbers
   * \param c New coloumn numbers
   */
  void resize(int r, int c) noexcept {
    if (r * c != size()) {
      if (r * c > _capacity) {
        if (dataPtr != nullptr) {
          if constexpr (std::is_trivially_destructible_v<Scalar_t>)
            for (int i = 0; i < _capacity; i++) {
              // alloc().destroy(dataPtr + i);
              (dataPtr + i)->~Scalar_t();
            }
          alloc().deallocate(dataPtr, _capacity);
        }

        dataPtr = alloc().allocate(r * c);
        _capacity = r * c;

        if constexpr (std::is_trivially_default_constructible_v<Scalar_t>)
          for (int i = 0; i < _capacity; i++) {
            // alloc().construct(dataPtr + i);
            dataPtr[i] = Scalar_t();
          }
      }
    }

    rowNum = r;
    colNum = c;
  }

  /**
   * \brief Same API to eigen's matrices
   */
  inline int rows() const noexcept { return rowNum; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline int cols() const noexcept { return colNum; }

  /**
   * \brief Same API to std::vector
   *
   */
  inline int capacity() const noexcept { return _capacity; }

  /**
   * \brief Same API to std::vector
   */
  inline const Scalar_t &operator[](int n) const noexcept { return dataPtr[n]; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline const Scalar_t &operator()(int n) const noexcept { return dataPtr[n]; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline const Scalar_t &operator()(int r, int c) const noexcept { return dataPtr[rowNum * c + r]; }

  /**
   * \brief Same API to std::vector
   */
  inline Scalar_t &operator[](int n) noexcept { return dataPtr[n]; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline Scalar_t &operator()(int n) noexcept { return dataPtr[n]; }

  /**
   * \brief Same API to eigen's matrices
   */
  inline Scalar_t &operator()(int r, int c) noexcept { return dataPtr[rowNum * c + r]; }

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

  // static constexpr bool isClass = std::is_class<Scalar_t>::value;
  static const bool isFixedSize = false;
  static constexpr int rowsAtCompileTime = Eigen::Dynamic;
  static constexpr int colsAtCompileTime = Eigen::Dynamic;

 protected:
  Scalar_t *dataPtr;
  int rowNum;
  int colNum;
  int _capacity;

 private:
  inline static allocator_t &alloc() noexcept {
    static allocator_t alloctor;
    return alloctor;
  }
};

}  // namespace heu

#endif  // MATRIXDYNAMICSIZE_H
