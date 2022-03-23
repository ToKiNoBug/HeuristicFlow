// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_BOXREAL_H
#define EIGEN_HEU_BOXREAL_H

#include "InternalHeaderCheck.h"
#include "BoxRTCTDims.hpp"

namespace Eigen {

namespace internal {

/**
 * \ingroup HEU_EAGlobal
 * \class RealBoxBase
 * \brief Internal base class for all kinds of real boxes.
 *
 * \tparam Scalar_t Type of scalar
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container
 * \tparam BS Boxshape
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Min value (DivCode)
 * \tparam MaxCT Max value (DivCode)
 *
 * \note Real boxes of different types
 */
template <typename Scalar_t, int Dim, DoubleVectorOption DVO, BoxShape BS, bool isFixedRange,
          TemplateVal_t<Scalar_t> MinCT, TemplateVal_t<Scalar_t> MaxCT>
class RealBoxBase : public BoxDims<Scalar_t, Dim, DVO, BS, isFixedRange, MinCT, MaxCT> {
 private:
  static_assert(std::is_floating_point<Scalar_t>::value, "Scalar_t must be a floating point number");
  static_assert(Dim > 0 || Dim == Eigen::Dynamic, "Invalid template parameter Dim");

 public:
  /**
   * \brief Type of encoding.
   *
   */
  static const constexpr EncodeType Encoding = EncodeType::Real;
};

/**
 * \ingroup HEU_EAGlobal
 * \class learnRateBody
 * \brief Internal base class that for runtime learning rate.
 *
 * \tparam Var_t Type of floating-point number or vector.
 */
template <typename Var_t>
struct learnRateBody {
 protected:
  Var_t _learnRate;  ///< Valueu of learn rate, a scalar if squarebox and vector for non-square box

 public:
  inline Var_t& learnRate() { return _learnRate; }  ///< Non-const reference to _learnRate

  inline const Var_t& learnRate() const { return _learnRate; }  ///< Const reference to learnRate

  inline void setLearnRate(const Var_t& v) { _learnRate = v; }  ///< Change the value of learn rate
};

/**
 * \ingroup HEU_EAGlobal
 * \class RealBox
 * \brief Compile-time ranged box with fixed learning rate
 *
 * \tparam Scalar_t Type of scalar
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container
 * \tparam BS Boxshape
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Min value (DivCode)
 * \tparam MaxCT Max value (DivCode)
 * \tparam LearnRateCT Learn rate at compile time(DivCode)
 */
template <typename Scalar_t, int Dim, DoubleVectorOption DVO, BoxShape BS, bool isFixedRange = false,
          TemplateVal_t<Scalar_t> MinCT = TemplateVal_t<Scalar_t>(1),
          TemplateVal_t<Scalar_t> MaxCT = TemplateVal_t<Scalar_t>(1),
          TemplateVal_t<Scalar_t> LearnRateCT = TemplateVal_t<Scalar_t>(1)>
class RealBox : public RealBoxBase<Scalar_t, Dim, DVO, BS, isFixedRange, MinCT, MaxCT> {
 private:
  static_assert(isFixedRange, "Wrong specialization of RealBox");
  static const constexpr Scalar_t learnRateCT = DivDecode<LearnRateCT>::real;

 public:
  inline constexpr Scalar_t learnRate() const { return learnRateCT; }  ///< Learn rate at compile time
};

/**
 * @brief runtime ranged box with learning rate
 */

/**
 * \ingroup HEU_EAGlobal
 * \brief Runtime-ranged box with runtime range and learning rate
 *
 * \tparam Scalar_t Type of scalar
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container
 * \tparam BS Boxshape
 * \tparam MinCT Meaningless for this specialization
 * \tparam MaxCT Meaningless for this specialization
 * \tparam LearnRateCT Meaningless for this specialization
 */
template <typename Scalar_t, int Dim, DoubleVectorOption DVO, BoxShape BS, TemplateVal_t<Scalar_t> MinCT,
          TemplateVal_t<Scalar_t> MaxCT, TemplateVal_t<Scalar_t> LearnRateCT>
class RealBox<Scalar_t, Dim, DVO, BS, false, MinCT, MaxCT, LearnRateCT>
    : public RealBoxBase<Scalar_t, Dim, DVO, BS, false, MinCT, MaxCT>,
      public learnRateBody<
          typename std::conditional<BS == BoxShape::RECTANGLE_BOX, Container<Scalar_t, Dim, DVO>, Scalar_t>::type> {};

}  // namespace internal

}  // namespace Eigen

#endif  // EIGEN_HEU_BOXREAL_H
