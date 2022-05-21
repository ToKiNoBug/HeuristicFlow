
#ifndef Heu_MATRIXMAP_HPP
#define Heu_MATRIXMAP_HPP

namespace heu {

#include "InternalHeaderCheck.h"

template <typename Scalar_t>
class MatrixMap {
 protected:
  using fast_t = typename std::conditional<sizeof(Scalar_t) <= 3 * sizeof(void *), Scalar_t, const Scalar_t &>::type;

 public:
  MatrixMap(Scalar_t *d, size_t _r, size_t _c) {
    p = d;
    r = _r;
    c = _c;
  }
  /**
   * @brief Construct a new Matrix Map object through shallow copying
   *
   * @param s
   */
  MatrixMap(const MatrixMap &s) {
    r = s.r;
    c = s.c;
    p = s.p;
  }

  ~MatrixMap() = default;

  void deepCopyTo(MatrixMap *dst) const {
    dst->r = r;
    dst->c = c;
    for (size_t i = 0; i < size(); i++) {
      dst->operator()(i) = operator()(i);
    }
  }

  inline size_t rows() const noexcept { return r; }

  inline size_t cols() const noexcept { return c; }

  inline size_t size() const noexcept { return r * c; }

  inline Scalar_t *data() noexcept { return p; }

  inline const Scalar_t *data() const noexcept { return p; }

  inline Scalar_t *begin() noexcept { return p; }

  inline Scalar_t *end() noexcept { return p + r * c; }

  inline const Scalar_t *cbegin() const noexcept { return p; }

  inline const Scalar_t *cend() const noexcept { return p + r * c; }

  inline void resize(size_t _r, size_t _c) const noexcept {
    r = _r;
    c = _c;
  }

  inline Scalar_t &operator()(size_t i) { return p[i]; }

  inline Scalar_t &operator[](size_t i) { return p[i]; }

  inline Scalar_t &operator()(size_t _r, size_t _c) { return p[r * _c + _r]; }

  inline const Scalar_t &operator()(size_t i) const { return p[i]; }

  inline const Scalar_t &operator[](size_t i) const { return p[i]; }

  inline const Scalar_t &operator()(size_t _r, size_t _c) const { return p[r * _c + _r]; }

  inline void fill(fast_t src) {
    for (auto &i : *this) {
      i = src;
    }
  }

  static const bool isFixedSize = false;

 protected:
  Scalar_t *p;
  size_t r;
  size_t c;
};

}  // namespace heu

#endif  //  Heu_MATRIXMAP_HPP