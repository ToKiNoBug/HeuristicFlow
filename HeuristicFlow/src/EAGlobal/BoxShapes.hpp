
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
