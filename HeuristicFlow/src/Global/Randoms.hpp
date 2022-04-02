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

#ifndef HEU_RANDOMS_HPP
#define HEU_RANDOMS_HPP

#include <stdint.h>
#include <random>
#include <cmath>
#include <chrono>

#include "InternalHeaderCheck.h"

namespace heu {

namespace internal {

#ifdef __GNUC__
#if (defined __WIN32) || (defined __WIN64)
//  MingW's implementation for std::random_device doesn't produce a real random number
#define HEU_std_random_device_NOT_RELIABLE
#endif
#endif  //#ifdef __CNUC__

/**
 * \ingroup HEU_GLOBAL
 * \brief Internal global std::random device
 *
 * \return std::random_device& A reference to this static variable
 */
inline std::random_device& global_random_device() {
  static std::random_device rdv;
  return rdv;
}
/**
 * \ingroup HEU_GLOBAL
 * \brief Internal global std::mt19937 used as a high-performance
 * random number generater.
 *
 * \return std::mt19937& A reference to this instance.
 */
inline std::mt19937& global_mt19937() {
#ifdef HEU_std_random_device_NOT_RELIABLE
  // Use a hash value of system clock as the random seed
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

/**
 * \ingroup HEU_GLOBAL
 * \brief Uniform random number (double) in range [0,1)
 *
 * \return double random number
 */
inline double randD() {
  static std::uniform_real_distribution<double> rnd(0, 1);
  return rnd(internal::global_mt19937());
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Uniform random number (double) in range [min,max)
 *
 * \param min Minimum value
 * \param max Maximum value
 * \return double random number
 */
inline double randD(const double min, const double max) { return (max - min) * randD() + min; }

/**
 * \ingroup HEU_GLOBAL
 * \brief Uniform random number (float) in range [0,1)
 *
 * \return double random number
 */
inline float randF() {
  static std::uniform_real_distribution<float> rnd(0, 1);
  return rnd(internal::global_mt19937());
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Uniform random index in range [0,size)
 *
 * \tparam int_t Type of integer
 * \param size Maxixmum size value
 * \return int_t Random index in range [0,size-1]
 */
template <typename int_t>
inline int_t randIdx(int_t size) {
  static_assert(std::is_integral<int_t>::value, "int_t must be integer");
  return int_t(randF() * size);
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Uniform random index in range [min,max_plus_1)
 *
 * \tparam int_t Type of integer
 * \param min Minimun value
 * \param max_plus_1 One plus the maximum value
 * \return int_t Random index in range [min,max_plus_1 -1]
 */
template <typename int_t>
inline int_t randIdx(int_t min, int_t max_plus_1) {
  static_assert(std::is_integral<int_t>::value, "int_t must be integer");
  return int_t((max_plus_1 - min) * randF() + min);
}

}  // namespace heu

#endif  // HEU_RANDOMS_HPP
