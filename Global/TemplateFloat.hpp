/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef OptimT_TEMPLATEFLOAT_HPP
#define OptimT_TEMPLATEFLOAT_HPP

#include <stdint.h>

namespace OptimT {

/**
 * @brief Encode double by division of int32 and uint32 stored in uint64
 * 
 */
enum DivCode : uint64_t {
    Half=4294967298,
    Pi=2088567404207137453,
    Pi_mul_2=993512999210214248,
    Pi_mul_3=2111843339415065010,
    Pi_mul_4=4223686678804044427,
    Pi_mul_6=2111843339388979417,
    Pi_div_2=1744352159520604142,
    Pi_div_3=458953659822443224,
    Pi_div_4=1055921669994474028,
    Pi_div_6=114738414981121388,
    Sqrt2=1464156825348615605,
    one_div_Sqrt2=566232695896337324,
    E=1518759938472606051,
};

/**
 * @brief Metafunction to encode the numerator and denominator into uint64
 * 
 * @tparam a numerator
 * @tparam b denominator
 */
template<int32_t a,uint32_t b>
struct encode
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
struct decode
{
public:
    constexpr static const int32_t numerator=int32_t(dc>>32);
    constexpr static const uint32_t denominator=dc&0xFFFFFFFF;
    constexpr static const double real=double(numerator)/denominator;
};


enum PowCode : uint64_t {

};



template<uint64_t val,uint64_t threshold=10000000000000000ULL>
struct amplifier
{
private:
static const constexpr bool need2amp=(val<threshold);
public:
static const constexpr uint64_t result=need2amp?(amplifier<10*val,threshold>::result):val;
};

template<uint64_t threshold>
struct amplifier<0,threshold>
{
static const constexpr uint64_t result=0;
};

template<int16_t pow>
struct OneE
{
private:
static const constexpr int16_t direction=(pow>0)?-1:+1;
public:
static const constexpr double result=
        ((pow>0)?10.0:0.1)
            *OneE<pow+direction>::result;
};

template<>
struct OneE<0>
{
static const constexpr double result=1;
};

//stores 16 diget(dec) of precision
//1 bit for sign,54 bits for significant ,1+8 bits for power,
template<int64_t significant,int16_t power>
struct powEncode
{
private:
static const constexpr uint64_t threshold=1e16;
static const constexpr bool isNegative=significant<0;
static const constexpr uint64_t absVal=(isNegative)?(-significant):significant;
static const constexpr bool isSigValid=(absVal/10<threshold);

static_assert(isSigValid,"Unsupported 16+ decimal digits for precision");


static const constexpr uint64_t recordedSignificant=
        amplifier<absVal,threshold>::result;

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
struct powDecode
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
    static const constexpr double powPart=OneE<pow>::result;
public:
    static const constexpr double real=digital*powPart;
};

}   //  namespace OptimT

#endif  //  OptimT_TEMPLATEFLOAT_HPP
