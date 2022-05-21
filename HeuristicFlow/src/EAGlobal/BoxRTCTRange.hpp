
#ifndef HEU_BOX_RTCTRange_HPP
#define HEU_BOX_RTCTRange_HPP

#include <stdint.h>
#include <assert.h>

#include "InternalHeaderCheck.h"
#include <HeuristicFlow/Global>

/**
 * This file has implementation for different min() and max()
 * and incomplete setMin() and setMax() implementation
 */

namespace heu {

namespace internal {

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class BoxDynamicRange
 * \brief Internal base class for various types of boxes.
 *
 * \tparam Scalar_t Type of scalars
 * \tparam Size Box dimensions
 * \tparam BS Boxshape, BoxShape::SQUARE_BOX or BoxShape::RECTANGLE_BOX
 * \tparam DVO Type of container.
 *
 * \note Non-square box with runtime range.
 */
template <typename Scalar_t, int Size, BoxShape BS, ContainerOption DVO>
class BoxDynamicRange {
 public:
  using Var_t = Container<Scalar_t, Size, DVO>;  ///< Type of decision variable
  static const constexpr bool isBox = true;      ///< Denotes this class to be a box-constriant
  static const constexpr char Flag[] = "Non-square box with runtime range.";  ///< Box type in string
  static const constexpr BoxShape Shape = BoxShape::RECTANGLE_BOX;            ///< Box type

  inline Var_t& min() { return _minV; }  ///< Get a non-const reference to minimum vector

  inline Var_t& max() { return _maxV; }  ///< Get a non-const reference to mmaximum vector

  inline const Var_t& min() const { return _minV; }  ///< Get a const reference to minimum vector

  inline const Var_t& max() const { return _maxV; }  ///< Get a const reference to minimum vector

  inline void setMin(const Var_t _m) { _minV = _m; }  ///< Set the minimum vector

  inline void setMax(const Var_t _m) { _maxV = _m; }  ///< Set the maximum vector

 protected:
  Var_t _minV;  ///< This member exists when it's non-square box
  Var_t _maxV;  ///< This member exists when it's non-square box

 private:
  static_assert(BS != BoxShape::SQUARE_BOX, "Non square box specialization used for square box");
};

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class BoxDynamicRange
 * \brief Internal base class for various types of boxes.
 *
 *
 * \tparam Scalar_t Type of scalars
 * \tparam Size Box dimensions
 * \tparam DVO Type of container.
 *
 * \note Square box with runtime range
 */
template <typename Scalar_t, int Size, ContainerOption DVO>
class BoxDynamicRange<Scalar_t, Size, BoxShape::SQUARE_BOX, DVO> {
 public:
  using Var_t = Container<Scalar_t, Size, DVO>;                           ///< Type of decision variable
  static const constexpr bool isBox = true;                               ///< Denotes this class to be a box-constriant
  static const constexpr char Flag[] = "Square box with runtime range.";  ///< Box type in string
  static const constexpr BoxShape Shape = BoxShape::SQUARE_BOX;           ///< Box type

  inline Scalar_t min() const { return _minS; }  ///< Minimum value for each dim

  inline Scalar_t max() const { return _maxS; }  ///< Maximum value for each dim

  inline Scalar_t& min() { return _minS; }  ///< Non-const reference to maximum value for each dim

  inline Scalar_t& max() { return _maxS; }  ///< Non-const reference to minimum value for each dim

  inline void setMin(Scalar_t s) { _minS = s; }  ///< Set maximum value for each dim

  inline void setMax(Scalar_t s) { _maxS = s; }  ///< Set minimum value for each dim

 protected:
  Scalar_t _minS;  ///< This member exists when it's a runtime square box
  Scalar_t _maxS;  ///< This member exists when it's a runtime square box
};

/**
 * \ingroup CXX14_METAHEURISTIC
 * \brief Types that can pass a Scalar_t in templates.
 *
 * \tparam Scalar_t Type of scalar
 * \return Sclar_t for integers, and DivCode for floating-point numbers.
 */
template <typename Scalar_t>
using TemplateVal_t = typename std::conditional<std::is_floating_point<Scalar_t>::value, DivCode, Scalar_t>::type;

/**
 * \ingroup CXX14_METAHEURISTIC
 * \class BoxFixedRange
 * \brief Internal base class for various types of boxes.
 *
 * \tparam Scalar_t Type of scalar
 * \tparam MinCT Minimum value at compile-time
 * \tparam MaxCT Maximum value at compile-time
 *
 * \note Square box with compile-time range
 */
template <typename Scalar_t, TemplateVal_t<Scalar_t> MinCT, TemplateVal_t<Scalar_t> MaxCT>
class BoxFixedRange {
 public:
  static const constexpr bool isBox = true;  ///< Denotes this class to be a box-constriant
  static const constexpr char Flag[] = "Square box with compile-time range";  ///< Box type in string
  static const constexpr BoxShape Shape = BoxShape::SQUARE_BOX;               ///< Box type

  /**
   * \brief Minimum value that is fixed at compile time.
   *
   * \return constexpr Scalar_t Minimum value
   */
  inline constexpr Scalar_t min() const { return Decoder<isFloatPoint>::minCT; }

  /**
   * \brief Maximum value that is fixed at compile time.
   *
   * \return constexpr Scalar_t Maximum value
   */
  inline constexpr Scalar_t max() const { return Decoder<isFloatPoint>::maxCT; }

 private:
  static const constexpr bool isFloatPoint = std::is_floating_point<Scalar_t>::value;

  template <bool isFl, typename unused = void>
  struct Decoder {
    static const constexpr Scalar_t minCT = MinCT;
    static const constexpr Scalar_t maxCT = MaxCT;
    static_assert(isFl == isFloatPoint, "Wrong specialization");
  };

  template <typename unused>
  struct Decoder<true, unused> {
    static const constexpr Scalar_t minCT = DivDecode<MinCT>::real;
    static const constexpr Scalar_t maxCT = DivDecode<MaxCT>::real;
    static_assert(isFloatPoint, "Wrong specialization");
  };
};

}  // namespace internal

}  // namespace heu

#endif  // HEU_BOX_RTCTRange_HPP
