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
#include <cmath>
#include <ctime>
#include <random>

namespace Heu {

///global random device(mt19937) for Heu
extern std::mt19937 global_mt19937;

///Calling anything in this namespace is deprecated
namespace HeuPrivate {    
inline uint32_t makeRandSeed() {
        static bool isFirstCalled=true;
        if(isFirstCalled) {
            uint32_t seed,seed_;
            if(sizeof(std::time_t)==4) {
                seed=uint32_t(std::time(nullptr));
            }
            else {
                uint64_t _s=std::time(nullptr);
                seed=(_s>>32)^(_s&0xFFFFFFFF);
            }
            seed_=seed;
            std::srand(seed);
            uint8_t * swapper=(uint8_t *)&seed;
            std::swap(swapper[0],swapper[3]);
            std::swap(swapper[1],swapper[2]);
            isFirstCalled=false;

            return seed_^seed;
        }
        else {
            return global_mt19937.operator()();
        }
    }
}   // HeuPrivate


///uniform random number in range [0,1)
inline double randD() {
    static std::uniform_real_distribution<double> rnd(0,1);
    return rnd(global_mt19937);
}

///uniform random number in range [min,max)
inline double randD(const double min,const double max) {
    return (max-min)*randD()+min;
}

}   // Heu

#endif // RANDOMS_HPP
