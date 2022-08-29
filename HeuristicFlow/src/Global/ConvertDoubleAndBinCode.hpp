
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

#ifndef HEU_CONVERTDOUBLEANDBINCODE_HPP
#define HEU_CONVERTDOUBLEANDBINCODE_HPP

#include "InternalHeaderCheck.h"

#include <cmath>
#include <limits>
#include <stdint.h>

#include "Enumerations.hpp"
namespace heu {
namespace {

template <typename float_t = double>
inline constexpr float_t expOf2AtCompileTime(int N) {
  float_t fl = 1.0;

  while (N > 11) {
    fl *= 1024;
    N -= 10;
  }

  while (N > 0) {
    fl *= 2;
    N--;
  }

  while (N < -11) {
    fl /= 1024;
    N += 10;
  }

  while (N < 0) {
    fl /= 2;
    N++;
  }

  return fl;
}
}  // namespace

/**
 * \ingroup HEU_GLOBAL
 * \brief The binary code to convey double through templates. This enum stores floatinig-point
 * numbers dierctly by their binary encoding (IEEE754 standard).
 *
 * \sa decode encode
 */
enum binCode64 : uint64_t {};

/**
 * \ingroup HEU_GLOBAL
 * \brief The binary code to convey float through templates. This enum stores floatinig-point
 * numbers dierctly by their binary encoding (IEEE754 standard).
 *
 * \sa decode encode
 */
enum binCode32 : uint32_t {};

/**
 * \ingroup HEU_GLOBAL
 * \brief Convert binary code to double
 *
 * \param code The binary code
 * \return constexpr double The floating-point number at compile time
 */
inline constexpr double decode(const binCode64 binCode) {
  const uint64_t code = uint64_t(binCode);
  constexpr uint64_t expMask = (0b11111111111ULL) << 52;

  constexpr uint64_t mantissaMask = (0b1ULL << 52) - 1ULL;

  if (((expMask & code) == expMask) && (mantissaMask & code) != 0) {
    //  if is nan
    return std::numeric_limits<double>::signaling_NaN();
  }

  const uint64_t absCode = code & ~(1ULL << 63);
  if (absCode >= (0b11111111111ULL << 52)) {  //  if is inf
    if (absCode == code) {
      return std::numeric_limits<double>::infinity();
    } else {
      return -std::numeric_limits<double>::infinity();
    }
  }

  double res = 0;
  const uint64_t exponent = (code & expMask) >> 52;

  {
    if (exponent != 0) res += 1;

    double delta = 0.5;
    uint64_t singleMantissaMask = 1ULL << 51;

    for (int i = 0; i < 52; i++) {
      if (code & singleMantissaMask) {
        res += delta;
      }

      delta /= 2;
      singleMantissaMask = singleMantissaMask >> 1;
    }
  }

  // return res;

  {
    if (exponent != 0) {
      int N = int(exponent) - 1023;
      while (N > 0) {
        res *= 2;
        N--;
      }

      while (N < 0) {
        res /= 2;
        N++;
      }
    } else {
      res /= expOf2AtCompileTime(1022);
    }

    // res *= expOf2AtCompileTime(int(exponent) - 1023);
  }

  // return res;

  const bool isNegative = code & (1ULL << 63);
  if (isNegative)
    return -res;
  else
    return res;
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Convert binary code to float
 *
 * \param code The binary code
 * \return constexpr double The floating-point number at compile time
 */
inline constexpr float decode(const binCode32 binCode) {
  const uint32_t code = uint32_t(binCode);
  constexpr uint32_t expMask = (0b11111111UL) << 23;

  constexpr uint32_t mantissaMask = (0b1UL << 23) - 1UL;

  if (((expMask & code) == expMask) && (mantissaMask & code) != 0) {
    //  if is nan
    return std::numeric_limits<float>::signaling_NaN();
  }

  const uint32_t absCode = code & ~(1UL << 31);
  if (absCode >= (0b11111111UL) << 23) {  //  if is inf
    if (absCode == code) {
      return std::numeric_limits<float>::infinity();
    } else {
      return -std::numeric_limits<float>::infinity();
    }
  }

  float res = 0;
  uint32_t exponent = (code & expMask) >> 23;

  {
    if (exponent != 0) res += 1;

    float delta = 0.5;
    uint32_t singleMantissaMask = 1UL << 22;

    for (int i = 0; i < 52; i++) {
      if (code & singleMantissaMask) {
        res += delta;
      }

      delta /= 2;
      singleMantissaMask = singleMantissaMask >> 1;
    }
  }

  // return res;

  {
    if (exponent != 0) {
      int N = int(exponent) - 127;
      while (N > 0) {
        res *= 2;
        N--;
      }

      while (N < 0) {
        res /= 2;
        N++;
      }
    } else {
      res /= (expOf2AtCompileTime<float>(126));
    }

    // res *= expOf2AtCompileTime(int(exponent) - 1023);
  }

  // return res;

  const bool isNegative = code & (1UL << 31);
  if (isNegative)
    return -res;
  else
    return res;
}

namespace {

constexpr bool isNotNegative(const double fl) { return fl >= 0; }

constexpr bool isNotNegative(const float fl) { return fl >= 0; }

constexpr uint64_t signCode(const double fl) { return uint64_t(!isNotNegative(fl)) << 63; }

constexpr uint32_t signCode(const float fl) { return uint32_t(!isNotNegative(fl)) << 31; }

constexpr double absVal(const double fl) { return (fl >= 0) ? (fl) : (-fl); }

constexpr float absVal(const float fl) { return (fl >= 0) ? (fl) : (-fl); }

constexpr uint64_t exponentVal(double fl) {
  uint64_t code = 0;
  fl = absVal(fl);

  constexpr double seper =
      decode(binCode64(0b0000000000001111111111111111111111111111111111111111111111111111));

  while (fl > 257 * seper) {
    fl /= 256;
    code += 8;
  }

  while (fl > seper) {
    fl /= 2;
    code++;
  }

  return code;
}

constexpr uint32_t exponentVal(float fl) {
  uint32_t code = 0;
  fl = absVal(fl);
  /*
constexpr double seper =
    decode(binCode64(0b0000000000001111111111111111111111111111111111111111111111111111));
    */
  constexpr float seper = decode(binCode32(0b00000000011111111111111111111111UL));

  while (fl > 257 * seper) {
    fl /= 256;
    code += 8;
  }

  while (fl > seper) {
    fl /= 2;
    code++;
  }

  return code;
}

constexpr uint64_t exponentCode(double fl) { return uint64_t(exponentVal(fl)) << 52; }

constexpr uint32_t exponentCode(float fl) { return uint32_t(exponentVal(fl)) << 23; }

constexpr double absValWithExponentCode0(double fl) {
  fl = absVal(fl);

  if (exponentVal(fl) == 0) {
    return fl;
  }

  return (fl * expOf2AtCompileTime(1023 - (int)exponentVal(fl)) - 1) * expOf2AtCompileTime(-1022);
}

constexpr float absValWithExponentCode0(float fl) {
  fl = absVal(fl);

  if (exponentVal(fl) == 0) {
    return fl;
  }

  return (fl * expOf2AtCompileTime<float>(127 - (int)exponentVal(fl)) - 1) *
         expOf2AtCompileTime<float>(-126);
}

constexpr uint64_t mantissaCode(double fl) {
  fl = absValWithExponentCode0(fl);
  uint64_t mask = 1ULL << 51;

  uint64_t result = 0;

  constexpr double compareNum = expOf2AtCompileTime(-1023);

  for (int i = 0; i < 52; i++) {
    if (fl >= compareNum) {
      result = result | mask;
      fl -= compareNum;
    }
    fl *= 2;
    mask = mask >> 1;
  }

  return result;
}

constexpr uint32_t mantissaCode(float fl) {
  fl = absValWithExponentCode0(fl);
  uint32_t mask = 1UL << 22;

  uint32_t result = 0;

  constexpr float compareNum = (expOf2AtCompileTime<float>(-127));

  for (int i = 0; i < 23; i++) {
    if (fl >= compareNum) {
      result = result | mask;
      fl -= compareNum;
    }
    fl *= 2;
    mask = mask >> 1;
  }

  return result;
}

}  // namespace

/**
 * \ingroup HEU_GLOBAL
 * \brief Encode a floaitng-point number to binCode64.
 *
 * \param d The floating point number
 * \return constexpr binCode64 The binary code.
 */
constexpr binCode64 encode(const double d) {
  if (

      !(d >= 0 || d <= 0)  //  if a number is neither negative nor non-negative, it is nan

  ) {
    return binCode64((0ULL) | (0b11111111111ULL << 52) | (1ULL << 51));
  }

  if (absVal(d) >= std::numeric_limits<double>::infinity()) {
    if (d > 0) {
      return binCode64(0b11111111111ULL << 52);
    } else {
      return binCode64(0b111111111111ULL << 52);
    }
  } else
    return binCode64(signCode(d) | exponentCode(d) | mantissaCode(d));
}

/**
 * \ingroup HEU_GLOBAL
 * \brief Encode a floaitng-point number to binCode32.
 *
 * \param d The floating point number
 * \return constexpr binCode32 The binary code.
 */
constexpr binCode32 encode(const float d) {
  if (

      !(d >= 0 || d <= 0)  //  if a number is neither negative nor non-negative, it is nan

  ) {
    return binCode32((0UL) | (0b11111111UL << 23) | (1UL << 22));
  }

  if (absVal(d) >= std::numeric_limits<float>::infinity()) {
    if (d > 0) {
      return binCode32(0b11111111111UL << 23);
    } else {
      return binCode32(0b111111111111UL << 23);
    }
  } else
    return binCode32(signCode(d) | exponentCode(d) | mantissaCode(d));
}

template <typename float_t>
using binCode_t = typename std::enable_if_t<
    std::is_floating_point_v<float_t>,
    std::conditional_t<std::is_same_v<float_t, float>, binCode32, binCode64>>;

}  // namespace heu

#endif  // HEU_CONVERTDOUBLEANDBINCODE_HPP