// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_TYPES_HPP
#define EIGEN_HEU_TYPES_HPP

#include <Eigen/Core>
#include <type_traits>
#include <array>
#include <vector>

#include "InternalHeaderCheck.h"
#include "Enumerations.hpp"
#include "Constants.hpp"

namespace Eigen {

/**
 * \ingroup HEU_Global
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
 * \ingroup HEU_Global
 * \brief C++ std container of double.
 *
 * \tparam Size Number of element. Use Eigen::Dynamic for dynamic size.
 */
template <int Size>
using stdVecD_t = stdContainer<double, Size>;

/**
 * \ingroup HEU_Global
 * \brief Container for different types, sizes and scalar types.
 *
 * \tparam scalar_t Type of elements.
 * \tparam Size Size at compile time. Use Eigen::Dynamic for dynamic size.
 * \tparam DVO Type of container.
 * Use DoubleVectorOption::Eigen for Eigen::Array<scalar_t,Size,1> and
 * other enum values for stdContainer<scalar_t,Size>.
 *
 */
template <typename scalar_t, int Size, DoubleVectorOption DVO>
using Container = typename std::conditional<(DVO != DoubleVectorOption::Eigen), stdContainer<scalar_t, Size>,
                                            Eigen::Array<double, Size, 1> >::type;

// template<DoubleVectorOption dvo,size_t Dim>
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

}  //   namespace Eigen

#endif  // EIGEN_HEU_TYPES_HPP
