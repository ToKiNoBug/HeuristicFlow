/*
 Copyright © 2021-2022  TokiNoBug
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

    if constexpr (std::is_trivially_copy_assignable_v<typename Derived::Scalar_t>) {
      memcpy(static_cast<Derived&>(*this).data(), src.data(),
             sizeof(typename Derived::Scalar_t) * src.size());
    } else {
      for (int idx = 0; idx < src.size(); idx++) {
        static_cast<Derived&>(*this).data()[idx] = src.data()[idx];
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