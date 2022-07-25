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

#ifndef HEU_TYPES_HPP
#define HEU_TYPES_HPP

#include <Eigen/Core>
#include <type_traits>
#include <array>
#include <vector>

#include "InternalHeaderCheck.h"
#include "Enumerations.hpp"
#include "Constants.hpp"

namespace heu {
/**
 * \ingroup HEU_GLOBAL
 * \brief C++ std container for fixed and dynamic sizes.

 * Use std::vector for dynamic size and std::array as fixed size.
 *
 * \tparam scalar_t Type of element
 * \tparam Dim Size at compile time. Use Eigen::Dynamic for dynamic size.
*/
template <typename scalar_t, int Dim>
using stdContainer = typename std::conditional<Dim == Eigen::Dynamic, std::vector<scalar_t>,
                                               std::array<scalar_t, size_t(Dim)> >::type;

#if __cplusplus >= 202002L
template <class T>
concept isVector = requires(T vec, int idx) {
  {vec[idx]};
  {vec.size()};
  {vec.data()};
};

template <class T>
concept EigenClass = requires(T mat, int r, int c, T matB) {
  isVector<T>;
  {mat(r)};
  {mat(r, c)};
  {mat.rows()};
  {mat.cols()};
  {mat.array()};
  {mat = matB};
  {mat += 1};
  {mat.setZero()};
  {mat.array() + matB.array()};
  {mat.array() - matB.array()};
  {mat.array() * matB.array()};
  {mat.array() / matB.array()};
  {mat.minCoeff()};
  {mat.maxCoeff()};
  {mat.minCoeff(&r)};
  {mat.maxCoeff(&r, &c)};
};

template <class Box>
concept isBoxConstraint = requires(Box b, int d, typename Box::Var_t v) {
  typename Box::Var_t;
  typename Box::Scalar_t;
  Box::Shape;
  isVector<typename Box::Var_t>;
  std::is_arithmetic_v<typename Box::Scalar_t>;
  {b.min()};
  {b.max()};
  {b.min(d)};
  {b.max(d)};
  {b.initialize(&v)};
  {b.initializeSize(&v)};
  {b.applyConstraint(&v)};
  {b.applyConstraint(&v, d)};
  {b.applyDelta(&v)};
  {b.dimensions()};
};

template <class Box>
concept isPSOBox = requires(Box b, int d, typename Box::Var_t v) {
  typename Box::Var_t;
  typename Box::Scalar_t;
  Box::Shape;
  isVector<typename Box::Var_t>;
  std::is_arithmetic_v<typename Box::Scalar_t>;
  std::is_floating_point_v<typename Box::Scalar_t>;
  {b.posMin()};
  {b.posMax()};
  {b.posMin(d)};
  {b.posMax(d)};
  {b.velocityMax()};
  {b.velocityMax(d)};
  {b.initialize(&v)};
  {b.initializeSize(&v)};
  {b.applyConstraint(&v)};
  {b.applyConstraint(&v, d)};
};

#endif  //  #if __cplusplus >=202002L

namespace internal {

template <int _ObjNum>
struct heu_initializeSize {
 public:
  template <typename A>
  inline static void resize(A* v, size_t sz) noexcept {
    assert(v->size() == sz);
  }
};

template <>
struct heu_initializeSize<Eigen::Dynamic> {
 public:
  template <typename A>
  inline static void resize(A* v, size_t sz) noexcept {
    v->resize(sz);
  }
};

}  //    namespace internal

namespace {

/*
 * \ingroup HEU_GLOBAL
 * \brief Get the type of element from a type of vector/matrix
 *
 * \tparam Vec_t Type of vector/matrix
 * \return type Type of element
 */
template <class Vec_t>
struct toElement {
  using type = typename std::decay<decltype(Vec_t()[0])>::type;
};

/*
 * \ingroup HEU_GLOBAL
 * \brief Determine whether T is a Eigen class
 *
 * \note A class that publicaly inherits from Eigen class is also taken as Eigen class, cause it
 * has Eigen APIs.
 *
 * \tparam T Type of vector/matrix
 * \return value Whether T is a Eigen class.
 */
template <typename T>
struct isEigenClass {
#if __cplusplus >= 202002L
  static constexpr bool value = ::heu::template EigenClass<T>;
#else
  static constexpr bool value =
      std::is_assignable<Eigen::ArrayXX<typename toElement<T>::type>, T>::value;
#endif  //#if __cplusplus >= 202002L
};

// get the compile-time size of a eigen class
template <class T>
struct getEigenClassSizeCT {
  static constexpr int value = T::SizeAtCompileTime;
  static constexpr int rowsCT = T::RowsAtCompileTime;
  static constexpr int colsCT = T::ColsAtCompileTime;
  static_assert(isEigenClass<T>::value, "T must be a Eigen class!");
};

// get the compile-time size of a std::array
template <class T>
struct stdArraySizeCT {
  static constexpr int value = sizeof(T) / sizeof(typename toElement<T>::type);
  static constexpr int rowsCT = value;
  static constexpr int colsCT = 1;
};

// get the compile-time size of a std::array / std::vector
template <class T>
struct getStdVectorOrArraySizeCT {
 private:
  using scalar_t = typename toElement<T>::type;
  static constexpr bool isSizeFixed = !std::is_same<std::vector<scalar_t>, T>::value;

 public:
  static constexpr int value = (isSizeFixed) ? (stdArraySizeCT<T>::value) : (Eigen::Dynamic);
  static constexpr int rowsCT = value;
  static constexpr int colsCT = 1;
};

template <class T, bool isTEigenClass>
struct getSizeCTOfAnyVector_v {
  // here implements when isTEigenClass==true
  static constexpr int value = getEigenClassSizeCT<T>::value;
  static constexpr int rowsCT = getEigenClassSizeCT<T>::rowsCT;
  static constexpr int colsCT = getEigenClassSizeCT<T>::colsCT;
};

template <class T>
struct getSizeCTOfAnyVector_v<T, false> {
  static constexpr int value = getStdVectorOrArraySizeCT<T>::value;
  static constexpr int rowsCT = getStdVectorOrArraySizeCT<T>::rowsCT;
  static constexpr int colsCT = getStdVectorOrArraySizeCT<T>::colsCT;
};

// get the compile-time size of a std vector /std array / eigen class
template <class T>
struct getSizeCTOfAnyVector {
  static constexpr bool isTEigenClass = isEigenClass<T>::value;
  static constexpr int value = getSizeCTOfAnyVector_v<T, isTEigenClass>::value;
  static constexpr int rowsCT = getSizeCTOfAnyVector_v<T, isTEigenClass>::rowsCT;
  static constexpr int colsCT = getSizeCTOfAnyVector_v<T, isTEigenClass>::colsCT;
};
}  // namespace

template <class T>
struct array_traits {
  using Scalar_t = typename toElement<T>::type;
#if __cplusplus >= 202002L
  static constexpr bool isEigenClass = ::heu::template EigenClass<T>;
#else
  static constexpr bool isEigenClass = ::heu::template isEigenClass<T>::value;
#endif  //#if __cplusplus >= 202002L
  static constexpr ContainerOption containerType =
      (isEigenClass) ? (ContainerOption::Eigen) : (ContainerOption::Std);
  static constexpr int sizeCT = getSizeCTOfAnyVector<T>::value;
  static constexpr int rowsCT = getSizeCTOfAnyVector<T>::rowsCT;
  static constexpr int colsCT = getSizeCTOfAnyVector<T>::colsCT;
  static constexpr bool isFixedSize = (sizeCT != Eigen::Dynamic);

  static constexpr bool isRowVector = rowsCT == 1;
  static constexpr bool isColVector = colsCT == 1;
  static constexpr bool isVector = (isRowVector) || (isColVector);

  static constexpr bool isMatrix = (!isRowVector) && (!isColVector);
};

template <class T, int size>
struct sizeMayMatch {
  static constexpr bool value =
      (!array_traits<T>::isFixedSize) || (array_traits<T>::sizeCT == size);
};

template <class Var_t>
typename array_traits<Var_t>::Scalar_t& at(Var_t& v, const int idx) {
  if constexpr (array_traits<Var_t>::isEigenClass) {
    return v(idx);
  } else {
    return v[idx];
  }
}

template <class Var_t>
typename array_traits<Var_t>::Scalar_t at(const Var_t& v, const int idx) {
  if constexpr (array_traits<Var_t>::isEigenClass) {
    return v(idx);
  } else {
    return v[idx];
  }
}

}  //   namespace heu

#endif  // HEU_TYPES_HPP
