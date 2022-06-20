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

#ifndef HEU_BOXINTERFACE_HPP
#define HEU_BOXINTERFACE_HPP

#include "InternalHeaderCheck.h"
#include "BoxRTCTDims.hpp"
#include "BoxReal.h"

namespace heu {

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Fixed(N) dim real(d,double) square(S) box
 *
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minval (DivCode), default value is 0
 * \tparam MaxCT Maxval (DivCode), default value is 0
 * \tparam LRCT Learn rate(DivCode), default value is 0
 *
 * Although Boxes has varities of implementation, they are common rules of their APIs:
 * 1. APIs that always exists:
 *  - `box.min() const` for minimum value (returns scalar for square box, vector otherwise)
 *  - `box.max() const` for maximum value (returns scalar for square box, vector otherwise)
 *  - `int box.dimensions() const` for box dimension number.
 *  - `typename BoxType::Var_t` is the type of decision variable.
 *
 * 2. APIs that only exists when box is not square box, or a square box which doesn't have a fixed range:
 *  - `box.setMin(min)` for setting minimum value (min should be a scalar for square box, and vector otherwise)
 *  - `box.setMax(max)` for setting maximum value (max should be a scalar for square box, and vector otherwise)
 *  - `box.min()` returns a non-const reference of minimum (scalar/vector)
 *  - `box.max()` returns a non-const reference of maximum (scalar/vector)
 *
 * 3. APIs that only exists when box's dimensions isn't fixed at compile time:
 *  - `box.setDimensions(int dim)` Set the box's dimension to dim.
 *
 * 4. APIs that only exists for real boxes:
 *  - `box.learnRate() const` returns a the learning rate of the real box (the word learning rate is borrowed from
 * neural network tranning, here it refers to the greatest mangnitude to change the value of a decision variable)
 *
 * 5. APIs that only exists for real boxes which isn't a square box or doesn't have a fixed range:
 *  - `box.setLearnRate(LR)` to set the value of LR.
 *  - `box.learnRate()` returns a non-const reference to the learning rate (scalar/vector)
 *
 */
template <int Dim, ContainerOption DVO = ContainerOption::Std, bool isFixedRange = false,
          internal::TemplateVal_t<double> MinCT = DivEncode<0, 1>::code,
          internal::TemplateVal_t<double> MaxCT = DivEncode<0, 1>::code,
          internal::TemplateVal_t<double> LRCT = DivEncode<0, 1>::code>
using BoxNdS = internal::RealBox<double, Dim, DVO, BoxShape::SQUARE_BOX, isFixedRange, MinCT, MaxCT, LRCT>;

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief runtime(X) dim real(d-double) square(S) box
 *
 * \tparam DVO Type of container.
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minval (DivCode), default value is 0
 * \tparam MaxCT Maxval (DivCode), default value is 0
 * \tparam LRCT Learn rate(DivCode), default value is 0
 */
template <ContainerOption DVO = ContainerOption::Std, bool isFixedRange = false,
          internal::TemplateVal_t<double> MinCT = DivEncode<0, 1>::code,
          internal::TemplateVal_t<double> MaxCT = DivEncode<0, 1>::code,
          internal::TemplateVal_t<double> LRCT = DivEncode<0, 1>::code>
using BoxXdS = internal::RealBox<double, Eigen::Dynamic, DVO, BoxShape::SQUARE_BOX, isFixedRange, MinCT, MaxCT, LRCT>;

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Fixed(N) dim real(d-double) nonsquare(N) box
 *
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 */
template <int Dim, ContainerOption DVO = ContainerOption::Std>
using BoxNdN = internal::RealBox<double, Dim, DVO, BoxShape::RECTANGLE_BOX>;

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief runtime(X) dim real(d-double) nonsquare(N) box
 *
 * \tparam DVO Type of container, std containers for defult value.
 */
template <ContainerOption DVO = ContainerOption::Std>
using BoxXdN = internal::RealBox<double, Eigen::Dynamic, DVO, BoxShape::RECTANGLE_BOX>;

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class BooleanBox
 * \brief Binary box
 *
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 */
template <int Dim, ContainerOption DVO = ContainerOption::Std>
class BooleanBox : public internal::BoxDims<bool, Dim, DVO, BoxShape::SQUARE_BOX, true, 0, 1> {
 private:
  static_assert(Dim > 0 || Dim == Eigen::Dynamic, "Invalid template parameter Dim");

 public:
  /**
   * \brief This member denotes that this class is of binary encoding
   *
   */
  static const constexpr EncodeType Encoding = EncodeType::Binary;
};

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Binary(b) box with fixed dims(N)
 *
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 */
template <int Dim, ContainerOption DVO = ContainerOption::Std>
using BoxNb = typename std::enable_if<Dim != Eigen::Dynamic, BooleanBox<Dim, DVO>>::type;

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Binary(b) box with dynamic(X) dims
 *
 * \tparam DVO Type of container, std containers as default.
 */
template <ContainerOption DVO = ContainerOption::Std>
using BoxXb = BooleanBox<Eigen::Dynamic, DVO>;

namespace internal {

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class SymbolBox
 * \brief Internal class for symbolic box
 *
 * \tparam Scalar_t Type of symbols
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 * \tparam BS Boxshape
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time
 * \tparam MaxCT Maximum value at compile time
 */
template <typename Scalar_t, int Dim, ContainerOption DVO = ContainerOption::Std, BoxShape BS = BoxShape::SQUARE_BOX,
          bool isFixedRange = false, Scalar_t MinCT = 0, Scalar_t MaxCT = 1>
class SymbolBox : public internal::BoxDims<Scalar_t, Dim, DVO, BS, isFixedRange, MinCT, MaxCT> {
 private:
  static_assert(std::is_integral<Scalar_t>::value, "Symbol box requires integer Scalar_t");
  static_assert(Dim > 0 || Dim == Eigen::Dynamic, "Invalid template parameter Dim");

 public:
  /**
   * \brief This member denotes that this box is of symbolic encoding
   *
   */
  static const constexpr EncodeType Encoding = EncodeType::Symbolic;
};

}  // namespace internal

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Square(S) symbolic(s) box with fixed dim(N)
 *
 * \tparam Scalar_t Type of symbols
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time, default val 0
 * \tparam MaxCT Maximum value at compile time, default val 1
 *
 */
template <typename Scalar_t, int Dim, ContainerOption DVO = ContainerOption::Std, bool isFixedRange = false,
          Scalar_t MinCT = 0, Scalar_t MaxCT = 1>
using BoxNsS =
    typename std::enable_if<Dim != Eigen::Dynamic, internal::SymbolBox<Scalar_t, Dim, DVO, BoxShape::SQUARE_BOX,
                                                                       isFixedRange, MinCT, MaxCT>>::type;

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Square(S) symbolic(s) box with runtime(X) dim
 *
 * \tparam Scalar_t Type of symbols
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container, std containers as default.
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time, default val 0
 * \tparam MaxCT Maximum value at compile time, default val 1
 */
template <typename Scalar_t, ContainerOption DVO = ContainerOption::Std, bool isFixedRange = false, Scalar_t MinCT = 0,
          Scalar_t MaxCT = 1>
using BoxXsS = internal::SymbolBox<Scalar_t, Eigen::Dynamic, DVO, BoxShape::SQUARE_BOX, isFixedRange, MinCT, MaxCT>;

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Non-square(N) symbolic(s) box with fixed(N) dim
 *
 * \tparam Scalar_t Type of symbols
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container, std containers as default.
 */
template <typename Scalar_t, int Dim, ContainerOption DVO = ContainerOption::Std>
using BoxNsN = typename std::enable_if<Dim != Eigen::Dynamic,
                                       internal::SymbolBox<Scalar_t, Dim, DVO, BoxShape::RECTANGLE_BOX>>::type;

/**
 * @brief Non-square symbolic box with runtime dim
 */

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Non-square(N) symbolic(s) box with dynamic(X) dim
 *
 * \tparam Scalar_t Type of symboxs
 * \tparam DVO Type of container, std containers as default.
 */
template <typename Scalar_t, ContainerOption DVO = ContainerOption::Std>
using BoxXsN = internal::SymbolBox<Scalar_t, Eigen::Dynamic, DVO, BoxShape::RECTANGLE_BOX>;

}  // namespace heu

#endif  // HEU_BOXINTERFACE_HPP
