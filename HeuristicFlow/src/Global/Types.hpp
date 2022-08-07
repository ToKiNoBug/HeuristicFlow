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

template <size_t bitNum>
using uintX_t = std::enable_if_t<
    bitNum <= 64,
    std::conditional_t<bitNum <= 8, uint8_t,
                       std::conditional_t<bitNum <= 16, uint16_t,
                                          std::conditional_t<bitNum <= 32, uint32_t, uint64_t>>>>;

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
                                               std::array<scalar_t, size_t(Dim)>>::type;

namespace {

template <class Vec>
struct __impl_isVector {
 private:
  template <class V>
  static decltype(V()[int()], V().size(),

                  std::enable_if<std::is_pointer_v<decltype(V().data())>>(),

                  int())
  fun(const V&) {
    return 0;
  }

  static void fun(...) { return; }

 public:
  static constexpr bool value = std::is_same_v<int, decltype(fun(Vec()))>;
};

}  // namespace

/**
 * \ingroup HEU_GLOBAL
 * \brief
 *
 * \tparam vec
 */
template <class vec>
constexpr bool isVector_v = __impl_isVector<vec>::value;

namespace {

template <class box>
struct __impl_isBoxConstraint {
 private:
  template <class B>
  static decltype(

      std::enable_if<std::is_arithmetic_v<typename B::Scalar_t>>(),

      std::enable_if<isVector_v<typename B::Var_t>>(),

      B::Shape,

      B().min(int()), B().max(int()),

      B().initialize((typename B::Var_t*)(nullptr)),
      B().applyConstraint((typename B::Var_t*)(nullptr)),
      B().applyConstraint((typename B::Var_t*)(nullptr), int()),
      B().applyDelta((typename B::Var_t*)(nullptr)),

      B().dimensions(),

      int())
  fun(const B&) {
    return 0;
  }

  static void fun(...) {}

 public:
  static constexpr bool value = std::is_same_v<int, decltype(fun(box()))>;
};

}  // namespace

template <class box>
constexpr bool isBoxConstraint_v = __impl_isBoxConstraint<box>::value;

namespace {

template <class Box>
struct __impl_isContinousBox {
 private:
  template <class B>
  static decltype(std::enable_if<isBoxConstraint_v<B>>(),

                  std::enable_if<std::is_floating_point_v<typename B::Scalar_t>>(),

                  B().delta(int()),

                  int())
  fun(const B&) {
    return 0;
  }

  static void fun(...) { return; }

 public:
  static constexpr bool value = std::is_same_v<int, decltype(fun(Box()))>;
};

}  // namespace

template <class Box>
constexpr bool isContinousBox_v = __impl_isContinousBox<Box>::value;

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
  {mat.array() += 1};
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
  {b.min(d)};
  {b.max(d)};
  {b.initialize(&v)};
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
 private:
  struct isEigen_t {};

  static bool fun(...) { return 0; }

  template <class U>
  static decltype(U()(0), U()(0, 0), U().rows(), U().cols(), U().array(), U().setZero(),
                  U().array() += U().array(), U().array() -= U().array(),
                  U().array() *= U().array(), U().array() /= U().array(), U().minCoeff(),
                  U().maxCoeff(), isEigen_t())
  fun(const U&) {
    return isEigen_t();
  }

 public:
#if __cplusplus >= 202002L
  static constexpr bool value = ::heu::template EigenClass<T>;
#else
  /*                                                                          \
     static constexpr bool value =                                                  \
                                                                                    \
         std::is_assignable<Eigen::ArrayXX<typename toElement<T>::type>, T>::value; \
         */
  static constexpr bool value = std::is_same_v<decltype(fun(T())), isEigen_t>;
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

  struct isNotStdArray_t {};

  template <class U>
  static decltype(U().resize(10), isNotStdArray_t()) fun(const U&) {
    return isNotStdArray_t();
  }

  static int fun(...) { return 0; }

 private:
  static constexpr bool isSizeFixed = !std::is_same_v<isNotStdArray_t, decltype(fun(T()))>;
  //! std::is_same<std::vector<scalar_t>, T>::value;

 public:
  static constexpr int value = (isSizeFixed) ? (stdArraySizeCT<T>::value) : (Eigen::Dynamic);
  static constexpr int rowsCT = value;
  static constexpr int colsCT = 1;
};

template <class T, bool isTEigenClass>
struct __impl_getSizeOfAnyVector {
  // here implements when isTEigenClass==true
  static constexpr int value = getEigenClassSizeCT<T>::value;
  static constexpr int rowsCT = getEigenClassSizeCT<T>::rowsCT;
  static constexpr int colsCT = getEigenClassSizeCT<T>::colsCT;
};

template <class T>
struct __impl_getSizeOfAnyVector<T, false> {
  static constexpr int value = getStdVectorOrArraySizeCT<T>::value;
  static constexpr int rowsCT = getStdVectorOrArraySizeCT<T>::rowsCT;
  static constexpr int colsCT = getStdVectorOrArraySizeCT<T>::colsCT;
};

// get the compile-time size of a std vector /std array / eigen class
template <class T>
struct getSizeCTOfAnyVector {
  static constexpr bool isTEigenClass = isEigenClass<T>::value;
  static constexpr int value = __impl_getSizeOfAnyVector<T, isTEigenClass>::value;
  static constexpr int rowsCT = __impl_getSizeOfAnyVector<T, isTEigenClass>::rowsCT;
  static constexpr int colsCT = __impl_getSizeOfAnyVector<T, isTEigenClass>::colsCT;
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
