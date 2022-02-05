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

}   //  namespace OptimT

#endif  //  OptimT_TEMPLATEFLOAT_HPP
