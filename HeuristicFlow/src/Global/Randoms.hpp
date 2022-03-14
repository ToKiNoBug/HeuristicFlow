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

#ifndef Heu_RANDOMS_HPP
#define Heu_RANDOMS_HPP

#include <stdint.h>
#include <random>
#include <cmath>
#include <ctime>

namespace Heu {

#ifdef __GNUC__
#if (defined __WIN32) || (defined __WIN64)
#define Heu_std_random_device_UNRELIABLE
#endif
#endif  //  #ifdef __CNUC__

inline std::mt19937 & global_mt19937();

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

///uniform random number in range [0,1)
inline double randD() {
    static std::uniform_real_distribution<double> rnd(0,1);
    return rnd(global_mt19937());
}

///uniform random number in range [min,max)
inline double randD(const double min,const double max) {
    return (max-min)*randD()+min;
}

}   // Heu

#endif // RANDOMS_HPP
