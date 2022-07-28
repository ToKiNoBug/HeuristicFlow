
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

inline constexpr double expOf2AtCompileTime(int N) {
  double fl = 1.0;

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

namespace {

constexpr bool isNotNegative(const double fl) { return fl >= 0; }

constexpr uint64_t signCode(const double fl) { return uint64_t(!isNotNegative(fl)) << 63; }

constexpr double absVal(const double fl) { return (fl >= 0) ? (fl) : (-fl); }

constexpr uint64_t exponentVal(double fl) {
  uint64_t code = 0;
  fl = absVal(fl);

  constexpr double seper =
      decode(binCode64(0b0000000000001111111111111111111111111111111111111111111111111111));
  while ((fl / seper) > 257) {
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

constexpr double absValWithExponentCode0(double fl) {
  fl = absVal(fl);

  if (exponentVal(fl) == 0) {
    return fl;
  }

  return (fl * expOf2AtCompileTime(1023 - exponentVal(fl)) - 1) * expOf2AtCompileTime(-1022);
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

}  // namespace heu

#endif  // HEU_CONVERTDOUBLEANDBINCODE_HPP