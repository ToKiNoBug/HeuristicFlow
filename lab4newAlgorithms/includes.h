/*
 Copyright © 2022  TokiNoBug
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

#ifndef HeuristicFlow_LAB_INCLUDES_H
#define HeuristicFlow_LAB_INCLUDES_H

#include <Eigen/Dense>
#include <HeuristicFlow/Global>
#include <HeuristicFlow/Genetic>
#include <HeuristicFlow/SimpleMatrix>


enum DivCode : uint64_t {
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

template<int32_t a,uint32_t b>
struct encode
{
public:
    constexpr static const uint64_t value=(uint64_t(a)<<32)|b;
    constexpr static const DivCode code=(DivCode)value;
};

inline DivCode encodeFun(int32_t a,uint32_t b) {
    return (DivCode)((uint64_t(a)<<32)|b);
}

template<int32_t a,uint32_t b>
struct division
{
public:
    constexpr static const double val=double(a)/b;
};

template<DivCode dc>
struct decode
{
    constexpr static const int32_t numerator=int32_t(dc>>32);
    constexpr static const uint32_t denominator=dc&0xFFFFFFFF;
    constexpr static const double real=division<numerator,denominator>::val;
};


template<uint64_t a,uint64_t b>
double power(double base) {
    div_t dt;
    return std::pow(base,division<a,b>::val);
}

template<DivCode dc>
double power(double base) {
    return std::pow(base,decode<dc>::real);
}

DivCode findCodeOf(const double d,double * error=nullptr);
#endif  //  HeuristicFlow_LAB_INCLUDES_H
