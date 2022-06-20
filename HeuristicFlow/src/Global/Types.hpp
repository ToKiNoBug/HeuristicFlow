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
 *
 * Use std::vector for dynamic size and std::array as fixed size.
 *
 * \tparam scalar_t Type of element
 * \tparam Dim Size at compile time. Use Eigen::Dynamic for dynamic size.
 */
template <typename scalar_t, int Dim>
using stdContainer =
    typename std::conditional<Dim == Eigen::Dynamic, std::vector<scalar_t>, std::array<scalar_t, Dim> >::type;

/**
 * \ingroup HEU_GLOBAL
 * \brief C++ std container of double.
 *
 * \tparam Size Number of element. Use Eigen::Dynamic for dynamic size.
 */
template <int Size>
using stdVecD_t = stdContainer<double, Size>;

/**
 * \ingroup HEU_GLOBAL
 * \brief Container for different types, sizes and scalar types.
 *
 * \tparam scalar_t Type of elements.
 * \tparam Size Size at compile time. Use Eigen::Dynamic for dynamic size.
 * \tparam DVO Type of container.
 * Use ContainerOption::Eigen for Eigen::Array<scalar_t,Size,1> and
 * other enum values for stdContainer<scalar_t,Size>.
 *
 */
template <typename scalar_t, int Size, ContainerOption DVO>
using Container = typename std::conditional<(DVO != ContainerOption::Eigen), stdContainer<scalar_t, Size>,
                                            Eigen::Array<scalar_t, Size, 1> >::type;

// template<ContainerOption dvo,size_t Dim>
// using FitnessVec_t= Container<double,Dim,dvo>;

namespace internal {

template <int _ObjNum>
struct heu_initializeSize {
 public:
  template <typename A>
  inline static void resize(A* v, size_t sz) {
    assert(v->size() == sz);
  }
};

template <>
struct heu_initializeSize<Eigen::Dynamic> {
 public:
  template <typename A>
  inline static void resize(A* v, size_t sz) {
    v->resize(sz);
  }
};

}  //    namespace internal

/**
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

/**
 * \ingroup HEU_GLOBAL
 * \brief Determine whether T is a Eigen class
 *
 * \note A class that publicaly inherits from Eigen class is also taken as Eigen class, cause it has Eigen APIs.
 *
 * \tparam T Type of vector/matrix
 * \return value Whether T is a Eigen class.
 */
template <typename T>
struct isEigenClass {
  static constexpr bool value = std::is_assignable<Eigen::ArrayXX<typename toElement<T>::type>, T>::value;
};

namespace internal {

// get the compile-time size of a eigen class
template <class T>
struct getEigenClassSizeCT {
  static constexpr int value = T::SizeAtCompileTime;
  static_assert(isEigenClass<T>::value, "T must be a Eigen class!");
};

// get the compile-time size of a std::array
template <class T>
struct stdArraySizeCT {
  static constexpr int value = sizeof(T) / sizeof(typename toElement<T>::type);
};

// get the compile-time size of a std::array / std::vector
template <class T>
struct getStdVectorOrArraySizeCT {
 private:
  using scalar_t = typename toElement<T>::type;
  static constexpr bool isSizeFixed = !std::is_same<std::vector<scalar_t>, T>::value;

 public:
  static constexpr int value = (isSizeFixed) ? (stdArraySizeCT<T>::value) : (Eigen::Dynamic);
};

template <class T, bool isTEigenClass>
struct getSizeCTOfAnyVector_v {
  // here implements when isTEigenClass==true
  static constexpr int value = getEigenClassSizeCT<T>::value;
};

template <class T>
struct getSizeCTOfAnyVector_v<T, false> {
  static constexpr int value = getStdVectorOrArraySizeCT<T>::value;
};

// get the compile-time size of a std vector /std array / eigen class
template <class T>
struct getSizeCTOfAnyVector {
  static constexpr bool isTEigenClass = isEigenClass<T>::value;
  static constexpr int value = getSizeCTOfAnyVector_v<T, isTEigenClass>::value;
};
}  // namespace internal

template <class T>
struct array_traits {
  using Scalar_t = typename toElement<T>::type;
  static constexpr bool isEigenClass = isEigenClass<T>::value;
  static constexpr ContainerOption containerType = (isEigenClass) ? (ContainerOption::Eigen) : (ContainerOption::Std);
  static constexpr int sizeCT = internal::getSizeCTOfAnyVector<T>::value;
};

}  //   namespace heu

#endif  // HEU_TYPES_HPP
