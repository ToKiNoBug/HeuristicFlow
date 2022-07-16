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

#ifndef HEU_BOXRTCTDIMS_HPP
#define HEU_BOXRTCTDIMS_HPP

#include "InternalHeaderCheck.h"
#include "BoxShapes.hpp"

namespace heu {

namespace internal {

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class BoxCTDim
 * \brief Internal base class for various types of boxes.
 *
 * \tparam Scalar_t Type of scalars
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container.
 * \tparam BS Boxshape
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time
 * \tparam MaxCT Maximum value at compile time
 *
 * \note Square box with fixed dims
 */
template <typename Scalar_t, int Dim, ContainerOption DVO, BoxShape BS, bool isFixedRange,
          TemplateVal_t<Scalar_t> MinCT, TemplateVal_t<Scalar_t> MaxCT>
class BoxCTDim : public SquareBox<Scalar_t, Dim, DVO, isFixedRange, MinCT, MaxCT> {
 private:
  static_assert(Dim != Eigen::Dynamic, "BoxFixedSize used for runtime sized box");
  static_assert(BS == BoxShape::SQUARE_BOX, "NonSquarebox used for squarebox");

  using Base_t = SquareBox<Scalar_t, Dim, DVO, isFixedRange, MinCT, MaxCT>;

 public:
  using Var_t = typename Base_t::Var_t;

  inline constexpr int dimensions() const noexcept { return Dim; }
};

/**
 * @brief NonSquarebox with fixed dims
 */
template <typename Scalar_t, int Dim, ContainerOption DVO, bool isFixedRange,
          TemplateVal_t<Scalar_t> MinCT, TemplateVal_t<Scalar_t> MaxCT>
class BoxCTDim<Scalar_t, Dim, DVO, BoxShape::RECTANGLE_BOX, isFixedRange, MinCT, MaxCT>
    : public NonsquareBox<Scalar_t, Dim, DVO> {
 private:
  static_assert(isFixedRange == false,
                "Compile-time box range for non-square box is not supported");
  static_assert(Dim != Eigen::Dynamic, "BoxFixedSize used for runtime sized box");
  using Base_t = NonsquareBox<Scalar_t, Dim, DVO>;

 public:
  using Var_t = typename Base_t::Var_t;

  inline constexpr int dimensions() const noexcept { return Dim; }
};

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class BoxRTDim
 * \brief Internal base class for various types of boxes.
 *
 * \tparam Scalar_t Type of scalar
 * \tparam DVO Type of container
 * \tparam BS Boxshape
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time
 * \tparam MaxCT Maximum value at compile time
 *
 * \note Square box with runtime dims
 */
template <typename Scalar_t, ContainerOption DVO, BoxShape BS, bool isFixedRange,
          TemplateVal_t<Scalar_t> MinCT, TemplateVal_t<Scalar_t> MaxCT>
class BoxRTDim : public SquareBox<Scalar_t, Eigen::Dynamic, DVO, isFixedRange, MinCT, MaxCT> {
 private:
  static_assert(BS == BoxShape::SQUARE_BOX, "NonSquarebox used for squarebox");
  using Base_t = SquareBox<Scalar_t, Eigen::Dynamic, DVO, isFixedRange, MinCT, MaxCT>;

 protected:
  int _dims;  ///< Exisits when it's a square box with dynamic dims

 public:
  BoxRTDim() { _dims = 0; }  ///< Set initial dimensions to be 0.

  /**
   * \brief Set the Dimensions object
   *
   * \param d Must be positive.
   */
  inline void setDimensions(int d) noexcept {
    assert(d > 0);
    _dims = d;
  }

  /**
   * \brief Get dimensions
   *
   * \return int The value of protected member _dims
   */
  inline int dimensions() const noexcept { return _dims; }
};

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class BoxRTDim
 * \brief Internal base class for various types of boxes.
 *
 * \tparam Scalar_t Type of scalar
 * \tparam DVO Type of container
 * \tparam isFixedRange  Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time, meanningless for this specialization
 * \tparam MaxCT Maximum value at compile time, meanningless for this specialization
 *
 * \note Non-square box with dynamic dims
 */
template <typename Scalar_t, ContainerOption DVO, bool isFixedRange, TemplateVal_t<Scalar_t> MinCT,
          TemplateVal_t<Scalar_t> MaxCT>
class BoxRTDim<Scalar_t, DVO, BoxShape::RECTANGLE_BOX, isFixedRange, MinCT, MaxCT>
    : public NonsquareBox<Scalar_t, Eigen::Dynamic, DVO> {
 private:
  static_assert(isFixedRange == false,
                "Compile-time box range for non-square box is not supported");
  using Base_t = NonsquareBox<Scalar_t, Eigen::Dynamic, DVO>;

 public:
  using Var_t = typename Base_t ::Var_t;

  /**
   * \brief Get dimensions
   *
   * \return int The size of _minV.
   * Runtime assertion will fail if size of _minV and _maxV differs.
   */
  inline int dimensions() const noexcept {
    assert(this->_minV.size() == this->_maxV.size());
    return this->_minV.size();
  }
};

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class BoxDims
 * \brief Internal conditional base class for boolean box and symbolic box.
 * It's also a base class of real box.
 *
 * \tparam Scalar_t Type of scalar
 * \tparam Dim Box dimensions
 * \tparam DVO Type of container
 * \tparam BS Boxshape
 * \tparam isFixedRange Whether the range is fixed at compile time
 * \tparam MinCT Minimum value at compile time
 * \tparam MaxCT Maximum value at compile time
 */
template <typename Scalar_t, int Dim, ContainerOption DVO, BoxShape BS, bool isFixedRange,
          TemplateVal_t<Scalar_t> MinCT, TemplateVal_t<Scalar_t> MaxCT>
class BoxDims
    : public std::conditional<Dim == Eigen::Dynamic,
                              BoxRTDim<Scalar_t, DVO, BS, isFixedRange, MinCT, MaxCT>,
                              BoxCTDim<Scalar_t, Dim, DVO, BS, isFixedRange, MinCT, MaxCT>>::type {
};

}  // namespace internal

}  // namespace heu

#endif  // HEU_BOXRTCTDIMS_HPP
