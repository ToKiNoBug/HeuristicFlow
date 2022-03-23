// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_BOXINTERFACE_HPP
#define EIGEN_HEU_BOXINTERFACE_HPP

#include "InternalHeaderCheck.h"
#include "BoxRTCTDims.hpp"
#include "BoxReal.h"

namespace Eigen {

/**
 * \ingroup HEU_EAGlobal
 * \brief Fixed(N) dim real(d,double) square(S) box
 *
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minval (DivCode), default value is 0
 * \tparam MaxCT Maxval (DivCode), default value is 0
 * \tparam LRCT Learn rate(DivCode), default value is 0
 */
template <int Dim, DoubleVectorOption DVO = DoubleVectorOption::Std, bool isFixedRange = false,
          internal::TemplateVal_t<double> MinCT = DivEncode<0, 1>::code,
          internal::TemplateVal_t<double> MaxCT = DivEncode<0, 1>::code,
          internal::TemplateVal_t<double> LRCT = DivEncode<0, 1>::code>
using BoxNdS = internal::RealBox<double, Dim, DVO, BoxShape::SQUARE_BOX, isFixedRange, MinCT, MaxCT, LRCT>;

/**
 *
 */

/**
 * \ingroup HEU_EAGlobal
 * \brief runtime(X) dim real(d-double) square(S) box
 *
 * \tparam DVO Type of container.
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minval (DivCode), default value is 0
 * \tparam MaxCT Maxval (DivCode), default value is 0
 * \tparam LRCT Learn rate(DivCode), default value is 0
 */
template <DoubleVectorOption DVO = DoubleVectorOption::Std, bool isFixedRange = false,
          internal::TemplateVal_t<double> MinCT = DivEncode<0, 1>::code,
          internal::TemplateVal_t<double> MaxCT = DivEncode<0, 1>::code,
          internal::TemplateVal_t<double> LRCT = DivEncode<0, 1>::code>
using BoxXdS = internal::RealBox<double, Eigen::Dynamic, DVO, BoxShape::SQUARE_BOX, isFixedRange, MinCT, MaxCT, LRCT>;

/**
 * \ingroup HEU_EAGlobal
 * \brief Fixed(N) dim real(d-double) nonsquare(N) box
 *
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 */
template <int Dim, DoubleVectorOption DVO = DoubleVectorOption::Std>
using BoxNdN = internal::RealBox<double, Dim, DVO, BoxShape::RECTANGLE_BOX>;

/**
 * \ingroup HEU_EAGlobal
 * \brief runtime(X) dim real(d-double) nonsquare(N) box
 *
 * \tparam DVO Type of container, std containers for defult value.
 */
template <DoubleVectorOption DVO = DoubleVectorOption::Std>
using BoxXdN = internal::RealBox<double, Eigen::Dynamic, DVO, BoxShape::RECTANGLE_BOX>;

/**
 * \ingroup HEU_EAGlobal
 * \class BooleanBox
 * \brief Binary box
 *
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 */
template <int Dim, DoubleVectorOption DVO = DoubleVectorOption::Std>
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
 * \ingroup HEU_EAGlobal
 * \brief Binary(b) box with fixed dims(N)
 *
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 */
template <int Dim, DoubleVectorOption DVO = DoubleVectorOption::Std>
using BoxNb = typename std::enable_if<Dim != Eigen::Dynamic, BooleanBox<Dim, DVO>>::type;

/**
 * \ingroup HEU_EAGlobal
 * \brief Binary(b) box with dynamic(X) dims
 *
 * \tparam DVO Type of container, std containers as default.
 */
template <DoubleVectorOption DVO = DoubleVectorOption::Std>
using BoxXb = BooleanBox<Eigen::Dynamic, DVO>;

namespace internal {

/**
 * \ingroup HEU_EAGlobal
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
template <typename Scalar_t, int Dim, DoubleVectorOption DVO = DoubleVectorOption::Std,
          BoxShape BS = BoxShape::SQUARE_BOX, bool isFixedRange = false, Scalar_t MinCT = 0, Scalar_t MaxCT = 1>
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
 * \ingroup HEU_EAGlobal
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
template <typename Scalar_t, int Dim, DoubleVectorOption DVO = DoubleVectorOption::Std, bool isFixedRange = false,
          Scalar_t MinCT = 0, Scalar_t MaxCT = 1>
using BoxNsS =
    typename std::enable_if<Dim != Eigen::Dynamic, internal::SymbolBox<Scalar_t, Dim, DVO, BoxShape::SQUARE_BOX,
                                                                       isFixedRange, MinCT, MaxCT>>::type;

/**
 * \ingroup HEU_EAGlobal
 * \brief Square(S) symbolic(s) box with runtime(X) dim
 *
 * \tparam Scalar_t Type of symbols
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container, std containers as default.
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time, default val 0
 * \tparam MaxCT Maximum value at compile time, default val 1
 */
template <typename Scalar_t, DoubleVectorOption DVO = DoubleVectorOption::Std, bool isFixedRange = false,
          Scalar_t MinCT = 0, Scalar_t MaxCT = 1>
using BoxXsS = internal::SymbolBox<Scalar_t, Eigen::Dynamic, DVO, BoxShape::SQUARE_BOX, isFixedRange, MinCT, MaxCT>;

/**
 * \ingroup HEU_EAGlobal
 * \brief Non-square(N) symbolic(s) box with fixed(N) dim
 *
 * \tparam Scalar_t Type of symbols
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container, std containers as default.
 */
template <typename Scalar_t, int Dim, DoubleVectorOption DVO = DoubleVectorOption::Std>
using BoxNsN = typename std::enable_if<Dim != Eigen::Dynamic,
                                       internal::SymbolBox<Scalar_t, Dim, DVO, BoxShape::RECTANGLE_BOX>>::type;

/**
 * @brief Non-square symbolic box with runtime dim
 */

/**
 * \ingroup HEU_EAGlobal
 * \brief Non-square(N) symbolic(s) box with dynamic(X) dim
 *
 * \tparam Scalar_t Type of symboxs
 * \tparam DVO Type of container, std containers as default.
 */
template <typename Scalar_t, DoubleVectorOption DVO = DoubleVectorOption::Std>
using BoxXsN = internal::SymbolBox<Scalar_t, Eigen::Dynamic, DVO, BoxShape::RECTANGLE_BOX>;

}  // namespace Eigen

#endif  // EIGEN_HEU_BOXINTERFACE_HPP
