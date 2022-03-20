// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_TEMPLATEFLOAT_HPP
#define Heu_TEMPLATEFLOAT_HPP

#include <stdint.h>

namespace Eigen
{

/**
 * @brief Encode double by division of int32 and uint32 stored in uint64
 * 
 */
enum DivCode : uint64_t {
    DivCode_Half=4294967298,
    DivCode_Pi=2088567404207137453,
    DivCode_Pi_mul_2=993512999210214248,
    DivCode_Pi_mul_3=2111843339415065010,
    DivCode_Pi_mul_4=4223686678804044427,
    DivCode_Pi_mul_6=2111843339388979417,
    DivCode_Pi_div_2=1744352159520604142,
    DivCode_Pi_div_3=458953659822443224,
    DivCode_Pi_div_4=1055921669994474028,
    DivCode_Pi_div_6=114738414981121388,
    DivCode_Sqrt2=1464156825348615605,
    DivCode_one_div_Sqrt2=566232695896337324,
    DivCode_E=1518759938472606051,
};

/**
 * @brief Metafunction to encode the numerator and denominator into uint64
 * 
 * @tparam a numerator
 * @tparam b denominator
 */
template<int32_t a,uint32_t b>
struct DivEncode
{
private:   
    constexpr static const uint64_t value=(uint64_t(a)<<32)|b;
public:
    constexpr static const DivCode code=(DivCode)value;
};

/**
 * @brief Metafunction to decode a DivCode back to the numerator and denominator and corresponding floating-point number
 * 
 * @tparam dc DivCode waiting to be unpacked.
 */
template<DivCode dc>
struct DivDecode
{
public:
    constexpr static const int32_t numerator=int32_t(dc>>32);
    constexpr static const uint32_t denominator=dc&0xFFFFFFFF;
    constexpr static const double real=double(numerator)/denominator;
};


enum PowCode : uint64_t {

};

namespace internal
{

template<uint64_t val,uint64_t threshold=10000000000000000ULL>
struct PowEncode_amplifier
{
private:
static const constexpr bool need2amp=(val<threshold);
public:
static const constexpr uint64_t result=need2amp?(PowEncode_amplifier<10*val,threshold>::result):val;
};

template<uint64_t threshold>
struct PowEncode_amplifier<0,threshold>
{
static const constexpr uint64_t result=0;
};

template<int16_t pow>
struct PowEncode_OneE
{
private:
static const constexpr int16_t direction=(pow>0)?-1:+1;
public:
static const constexpr double result=
        ((pow>0)?10.0:0.1)
            *PowEncode_OneE<pow+direction>::result;
};

template<>
struct PowEncode_OneE<0>
{
static const constexpr double result=1;
};

}   //  internal

//stores 16 diget(dec) of precision
//1 bit for sign,54 bits for significant ,1+8 bits for power,
template<int64_t significant,int16_t power>
struct PowEncode
{
private:
static const constexpr uint64_t threshold=1e16;
static const constexpr bool isNegative=significant<0;
static const constexpr uint64_t absVal=(isNegative)?(-significant):significant;
static const constexpr bool isSigValid=(absVal/10<threshold);

static_assert(isSigValid,"Unsupported 16+ decimal digits for precision");


static const constexpr uint64_t recordedSignificant=
        internal::PowEncode_amplifier<absVal,threshold>::result;

static const constexpr bool isPowNegative=power<0;
static_assert (power<=255,"Power shouldn't exceeds 255");
static_assert (power>=-255,"Power shouldn't be less than -255");
static const constexpr uint8_t absPowVal=isPowNegative?(-power):power;

static const constexpr uint64_t upperPart
    =(uint64_t(isNegative)<<63)|(recordedSignificant<<9);
static const constexpr uint64_t lowerPart
    =(uint64_t(isPowNegative)<<8)|(absPowVal);

public:
static const constexpr PowCode code=PowCode(upperPart|lowerPart);
};



template<PowCode pc>
struct PowDecode
{
private:
    static const uint64_t code=pc;
    static const uint64_t fullMask=~(0ULL);
    static const constexpr uint64_t upperSigMask=((fullMask<<1)>>10)<<9;
    static const constexpr uint64_t upperBitMask=1ULL<<63;
    static const constexpr bool isNegative=upperBitMask&code;
    static const constexpr int64_t absSig=(code&upperSigMask)>>9;
    static const constexpr int64_t sig=(isNegative?(-absSig):absSig);

    static const constexpr uint64_t lowerBitMask=1ULL<<8;
    static const constexpr uint64_t lowerPowMask=0xFF;
    static const constexpr bool isPowNegative=lowerBitMask&code;
    static const constexpr int16_t absPow=lowerPowMask&code;
    static const constexpr int16_t pow=(isPowNegative?(-absPow):absPow);
    static const constexpr double digital=(sig)/(1e16);
    static const constexpr double powPart=internal::PowEncode_OneE<pow>::result;
public:
    static const constexpr double real=digital*powPart;
};

}   //  namespace Eigen

#endif  //  Heu_TEMPLATEFLOAT_HPP
