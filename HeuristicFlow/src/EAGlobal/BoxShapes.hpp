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

#ifndef HEU_BOXSHAPES_HPP
#define HEU_BOXSHAPES_HPP

#include <Eigen/Dense>
#include "InternalHeaderCheck.h"
#include "BoxRTCTRange.hpp"

namespace heu {

namespace internal {

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class NonsquareBox
 * \brief Internal typedef of base class for various types of boxes.
 *
 * \tparam Scalar_t Type of scalar
 * \tparam Size Box dimensions
 * \tparam DVO Type of container
 *
 * \note Non-square box types.
 */
template <typename Scalar_t, int Size, ContainerOption DVO>
class NonsquareBox : public BoxDynamicRange<Scalar_t, Size, BoxShape::RECTANGLE_BOX, DVO> {};

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class SquareBox
 * \brief Internal base class for various types of boxes.
 *
 * \tparam Scalar_t Type of scalar
 * \tparam Size Box dimensions
 * \tparam DVO Type of container
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time
 * \tparam MaxCT Maximum value at compile time
 *
 * \note Square box with compile time range
 */
template <typename Scalar_t, int Size, ContainerOption DVO, bool isFixedRange, TemplateVal_t<Scalar_t> MinCT,
          TemplateVal_t<Scalar_t> MaxCT>
class SquareBox : public BoxFixedRange<Scalar_t, MinCT, MaxCT> {
 private:
  using Base_t = BoxFixedRange<Scalar_t, MinCT, MaxCT>;

 public:
  using Var_t = Container<Scalar_t, Size, DVO>;  ///< Type of decision variable
};

/**
 * @brief
 */

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class SquareBox
 * \brief Internal base class for various types of boxes.
 *
 * \tparam Scalar_t Type of scalar
 * \tparam Size Box dimensions
 * \tparam DVO Type of container
 * \tparam MinCT Minimum value at compile time, meanningless for this specialization
 * \tparam MaxCT Maximum value at compile time, meanningless for this specialization
 *
 * \note Square box with runtime ranges
 */
template <typename Scalar_t, int Size, ContainerOption DVO, TemplateVal_t<Scalar_t> MinCT,
          TemplateVal_t<Scalar_t> MaxCT>
class SquareBox<Scalar_t, Size, DVO, false, MinCT, MaxCT>
    : public BoxDynamicRange<Scalar_t, Size, BoxShape::SQUARE_BOX, DVO> {
 private:
  using Base_t = BoxDynamicRange<Scalar_t, Size, BoxShape::SQUARE_BOX, DVO>;

 public:
  using Var_t = typename Base_t::Var_t;
};

}  // namespace internal

}  // namespace heu

#endif  // HEU_BOXSHAPES_HPP
