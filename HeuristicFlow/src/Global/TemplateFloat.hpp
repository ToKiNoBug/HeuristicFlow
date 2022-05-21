

#ifndef HEU_TEMPLATEFLOAT_HPP
#define HEU_TEMPLATEFLOAT_HPP

#include <stdint.h>

#include "InternalHeaderCheck.h"

namespace heu {

/**
 * \ingroup HEU_GLOBAL
 * \brief This enumeration encode a floating-point number by a division of int32 and uint32 stored in uint64
 *
 * In this way, floating-point values can be passed through template parameters and knowen at compile time under C++20.
 *
 * Members of this enumeration are some fundamental math constants.
 *
 * \note This enum should be taken as integers instead of enumerations. The reason it's designed to be a enum is to
 * prevent from being mixed up with general integers -- when you misuse integers as DivCode without explicitly
 * converting its type, you receive an error.
 *
 * \sa PowCode DivEncode DivDecode
 *
 */
enum DivCode : uint64_t {
  DivCode_Half = 4294967298,
  DivCode_Pi = 2088567404207137453,
  DivCode_Pi_mul_2 = 993512999210214248,
  DivCode_Pi_mul_3 = 2111843339415065010,
  DivCode_Pi_mul_4 = 4223686678804044427,
  DivCode_Pi_mul_6 = 2111843339388979417,
  DivCode_Pi_div_2 = 1744352159520604142,
  DivCode_Pi_div_3 = 458953659822443224,
  DivCode_Pi_div_4 = 1055921669994474028,
  DivCode_Pi_div_6 = 114738414981121388,
  DivCode_Sqrt2 = 1464156825348615605,
  DivCode_one_div_Sqrt2 = 566232695896337324,
  DivCode_E = 1518759938472606051,
};

/**
 * \ingroup HEU_GLOBAL
 * \struct DivEncode
 * \brief Metafunction to encode the numerator and denominator into uint64.
 *
 *
 * In this way, floating-point values can be passed through template parameters
 * and knowen at compile time under C++20.
 *
 * \tparam a Numerator, it can be positive, negative or 0
 * \tparam b Denominator, not less than 1.
 * \return `DivCode` code the encoded
 *
 * \sa DivCode DivDecode
 */
template <int32_t a, uint32_t b>
struct DivEncode {
 private:
  constexpr static const uint64_t value = (uint64_t(a) << 32) | b;

 public:
  constexpr static const DivCode code = (DivCode)value;
};

/**
 * \ingroup HEU_GLOBAL
 * \struct DivDecode
 * \brief Metafunction to decode a DivCode back to the numerator and denominator
 * and corresponding floating-point number.
 *
 * \tparam `dc` `DivCode` to be decoded.
 * \return `double` real the decoded floating-point number at compile time.
 *
 * \sa DivCode DivEncode
 */
template <DivCode dc>
struct DivDecode {
 public:
  constexpr static const int32_t numerator = int32_t(dc >> 32);
  constexpr static const uint32_t denominator = dc & 0xFFFFFFFF;
  constexpr static const double real = double(numerator) / denominator;
};

/**
 * \ingroup HEU_GLOBAL
 * \brief This enumeration encode a floating-point number
 * by storing its sign, its coefficient and exponent into a uint64.
 *
 *
 * It has the same fuction like DivCode but can represent a larger range.
 *
 * PowEncode has 16 diget(decimal) of precision.
 *
 * In all the 64 binary bits, 1 for sign, 54 for significant digits and 9 bits for power.
 *
 * \sa DivCode PowEncode PowDecode
 */
enum PowCode : uint64_t {

};

//
// 1 bit for sign,54 bits for significant ,1+8 bits for power,

namespace internal {
// Multiplie val for 10 times until it gets greater than the threshold
template <uint64_t val, uint64_t threshold = 10000000000000000ULL>
struct PowEncode_amplifier {
 private:
  static const constexpr bool need2amp = (val < threshold);

 public:
  static const constexpr uint64_t result = need2amp ? (PowEncode_amplifier<10 * val, threshold>::result) : val;
};

template <uint64_t threshold>
struct PowEncode_amplifier<0, threshold> {
 public:
  static const constexpr uint64_t result = 0;
};

template <int16_t pow>
struct PowEncode_OneE {
 private:
  static const constexpr int16_t direction = (pow > 0) ? -1 : +1;

 public:
  static const constexpr double result = ((pow > 0) ? 10.0 : 0.1) * PowEncode_OneE<pow + direction>::result;
};

template <>
struct PowEncode_OneE<0> {
 public:
  static const constexpr double result = 1;
};

}  // namespace internal

/**
 * \ingroup HEU_GLOBAL
 * \struct PowEncode
 * \brief Metafunction to encode a floating-point number like scientific notation
 *
 * \tparam significant The coefficient without the dot.
 * \tparam power The exponet power.
 *
 * \return PowCode code The encoded powcode.
 *
 * \code {.cpp}
 * PowEncode<1919810,-20>::code==PowEncode<191981,-20>::code;  // 1.191981^-20.
 * \endcode
 *
 * \sa PowCode PowDecode
 *
 */
template <int64_t significant, int16_t power>
struct PowEncode {
 private:
  static const constexpr uint64_t threshold = 1e16;
  static const constexpr bool isNegative = significant < 0;
  static const constexpr uint64_t absVal = (isNegative) ? (-significant) : significant;
  static const constexpr bool isSigValid = (absVal / 10 < threshold);

  static_assert(isSigValid, "Unsupported 16+ decimal digits for precision");

  static const constexpr uint64_t recordedSignificant = internal::PowEncode_amplifier<absVal, threshold>::result;

  static const constexpr bool isPowNegative = power < 0;
  static_assert(power <= 255, "Power shouldn't exceeds 255");
  static_assert(power >= -255, "Power shouldn't be less than -255");
  static const constexpr uint8_t absPowVal = isPowNegative ? (-power) : power;

  static const constexpr uint64_t upperPart = (uint64_t(isNegative) << 63) | (recordedSignificant << 9);
  static const constexpr uint64_t lowerPart = (uint64_t(isPowNegative) << 8) | (absPowVal);

 public:
  static const constexpr PowCode code = PowCode(upperPart | lowerPart);
};

/**
 * \ingroup HEU_GLOBAL
 * \struct PowDecode
 * \brief Metafunction to decode a divcode back into its real numer.
 *
 * \tparam pc Powcode to be decoded.
 * \return double real The decoded floating-point number.
 *
 * \code {.cpp}
 * constexpr PowCode pCode=PowEncode<-123456,16>::code;
 * constexpr double real=PowDecode<pCode>::real;
 * \endcode
 *
 * \sa PowCode PowEncode
 */
template <PowCode pc>
struct PowDecode {
 private:
  static const uint64_t code = pc;
  static const uint64_t fullMask = ~(0ULL);
  static const constexpr uint64_t upperSigMask = ((fullMask << 1) >> 10) << 9;
  static const constexpr uint64_t upperBitMask = 1ULL << 63;
  static const constexpr bool isNegative = upperBitMask & code;
  static const constexpr int64_t absSig = (code & upperSigMask) >> 9;
  static const constexpr int64_t sig = (isNegative ? (-absSig) : absSig);

  static const constexpr uint64_t lowerBitMask = 1ULL << 8;
  static const constexpr uint64_t lowerPowMask = 0xFF;
  static const constexpr bool isPowNegative = lowerBitMask & code;
  static const constexpr int16_t absPow = lowerPowMask & code;
  static const constexpr int16_t pow = (isPowNegative ? (-absPow) : absPow);
  static const constexpr double digital = (sig) / (1e16);
  static const constexpr double powPart = internal::PowEncode_OneE<pow>::result;

 public:
  static const constexpr double real = digital * powPart;
};

}  // namespace heu

#endif  // HEU_TEMPLATEFLOAT_HPP
