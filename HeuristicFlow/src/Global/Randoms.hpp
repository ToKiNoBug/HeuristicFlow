// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_RANDOMS_HPP
#define EIGEN_HEU_RANDOMS_HPP

#include <stdint.h>
#include <random>
#include <cmath>
#include <chrono>

#include "InternalHeaderCheck.h"

namespace Eigen {

namespace internal {
/**
 * \brief MingW's implementation for std::random_device doesn't produce a real random number
 *
 */
#ifdef __GNUC__
#if (defined __WIN32) || (defined __WIN64)
#define EIGEN_HEU_std_random_device_NOT_RELIABLE
#endif
#endif  //#ifdef __CNUC__

/**
 * \brief Internal global std::random device
 *
 * \return std::random_device& A reference to this static variable
 */
inline std::random_device& global_random_device() {
  static std::random_device rdv;
  return rdv;
}
/**
 * \brief Internal global std::mt19937
 *
 * \return std::mt19937& Used as a high-performance random number generater.
 */
inline std::mt19937& global_mt19937() {
#ifdef EIGEN_HEU_std_random_device_NOT_RELIABLE
  static auto now = std::chrono::system_clock::now();
  static std::time_t time = std::chrono::system_clock::to_time_t(now);
  static uint32_t seed = std::hash<std::time_t>()(time);
  static std::mt19937 mt(seed);
#else
  static std::mt19937 mt(global_random_device()());
#endif

  return mt;
}

}  // namespace internal

/// uniform random number in range [0,1)
inline double ei_randD() {
  static std::uniform_real_distribution<double> rnd(0, 1);
  return rnd(internal::global_mt19937());
}

/// uniform random number in range [min,max)
inline double ei_randD(const double min, const double max) { return (max - min) * ei_randD() + min; }

inline float ei_randF() {
  static std::uniform_real_distribution<float> rnd(0, 1);
  return rnd(internal::global_mt19937());
}

template <typename int_t>
inline int_t ei_randIdx(int_t size) {
  static_assert(std::is_integral<int_t>::value, "int_t must be integer");
  return int_t(ei_randF() * size);
}

template <typename int_t>
inline int_t ei_randIdx(int_t min, int_t max_plus_1) {
  static_assert(std::is_integral<int_t>::value, "int_t must be integer");
  return int_t((max_plus_1 - min) * ei_randF() + min);
}

}  // namespace Eigen

#endif  // EIGEN_HEU_RANDOMS_HPP
