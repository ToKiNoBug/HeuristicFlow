// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_RANDOMS_HPP
#define Heu_RANDOMS_HPP

#include <stdint.h>
#include <random>
#include <cmath>
#include <ctime>

namespace Eigen
{

namespace internal
{

#ifdef __GNUC__
#if (defined __WIN32) || (defined __WIN64)
#define Heu_std_random_device_UNRELIABLE
#endif
#endif  //  #ifdef __CNUC__


inline std::random_device & global_random_device() {
    static std::random_device rdv;
    return rdv;
}

inline std::mt19937 & global_mt19937() {
#ifdef Heu_std_random_device_UNRELIABLE
    static std::time_t time=std::time(nullptr);
    static uint32_t seed=std::hash<std::time_t>()(time);
    static std::mt19937 mt(seed);
#else
    static std::mt19937 mt(global_random_device()());
#endif
    return mt;
}

}   //  internal

///uniform random number in range [0,1)
inline double randD() {
    static std::uniform_real_distribution<double> rnd(0,1);
    return rnd(internal::global_mt19937());
}

///uniform random number in range [min,max)
inline double randD(const double min,const double max) {
    return (max-min)*randD()+min;
}

inline float randF() {
    static std::uniform_real_distribution<float> rnd(0,1);
    return rnd(internal::global_mt19937());
}

template<typename int_t>
inline int_t randIdx(int_t size) {
    static_assert (std::is_integral<int_t>::value,"int_t must be integer");
    return int_t(randF()*size);
}

template<typename int_t>
inline int_t randIdx(int_t min,int_t max_plus_1) {
    static_assert (std::is_integral<int_t>::value,"int_t must be integer");
    return int_t((max_plus_1-min)*randF()+min);
}


}   //  namespace Eigen

#endif // RANDOMS_HPP
