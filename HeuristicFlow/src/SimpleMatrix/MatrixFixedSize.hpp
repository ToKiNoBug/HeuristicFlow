
#ifndef Heu_MATRIXFIXEDSIZE_H
#define Heu_MATRIXFIXEDSIZE_H

#include <stdint.h>
#include <assert.h>
#include <array>

#include "InternalHeaderCheck.h"

namespace heu {

/**
 * \ingroup HEU_SIMPLEMATRIX
 * \brief Matrix with fixed size
 *
 * \tparam Scalar_t Type of element
 * \tparam Rows Row number
 * \tparam Cols Coloumn number
 */
template <class Scalar_t, size_t Rows, size_t Cols>
class MatrixFixedSize {
 protected:
  using fast_t = typename std::conditional<sizeof(Scalar_t) <= 3 * sizeof(void *), Scalar_t, const Scalar_t &>::type;

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
    for (size_t i = 0; i < size(); i++) {
      array[i] = src.array[i];
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
}  // namespace heu

#endif  // MATRIXFIXEDSIZE_H
