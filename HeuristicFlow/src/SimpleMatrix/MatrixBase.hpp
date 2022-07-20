#ifndef HEU_MATRIXBASE_HPP
#define HEU_MATRIXBASE_HPP

#include <type_traits>

#include "InternalHeaderCheck.h"

namespace heu {
template <class Derived>
class MatrixBase {
 protected:
  template <class DerivedB>
  Derived& operator=(const DerivedB& src) noexcept {
    static_assert(std::is_same_v<typename Derived::Scalar_t, typename DerivedB::Scalar_t>,
                  "You mixed different types of scalar.");
    // make size assertions and resize if this is dynamic-sized.
    if constexpr (Derived::isFixedSize) {
      if constexpr (DerivedB::isFixedSize) {
        static_assert(DerivedB::rowsAtCompileTime == Derived::rowsAtCompileTime &&
                          DerivedB::colsAtCompileTime == Derived::colsAtCompileTime,
                      "You mixed different sizes of matrices");
      } else {
        assert(Derived::rowsAtCompileTime == src.rows());
        assert(Derived::colsAtCompileTime == src.cols());
      }
    } else {
      static_cast<Derived&>(*this).resize(src.rows(), src.cols());
    }

    if constexpr (Derived::isClass) {
      memcpy(static_cast<Derived&>(*this).data(), src.data(),
             sizeof(typename Derived::Scalar_t) * src.size());
    } else {
      for (int idx = 0; idx < src.size(); idx++) {
        static_cast<Derived&>(*this).data()[idx] = src.data[idx];
      }
    }

    // static_cast<Derived&>(*this);
    return static_cast<Derived&>(*this);
  }
};
/*
template <class DerivedA, class DerivedB>
const auto& operator=(MatrixBase<DerviedA>& dst, const MatrixBase<DerivedB>& src) noexcept {
  ;
}
*/

}  // namespace heu

#endif  //  HEU_MATRIXBASE_HPP