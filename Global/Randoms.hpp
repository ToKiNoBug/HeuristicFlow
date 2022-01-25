#ifndef RANDOMS_HPP
#define RANDOMS_HPP

#include <stdint.h>
#include <cmath>
#include <ctime>
#include <random>

namespace OptimT {

///global random device(mt19937) for OptimT
extern std::mt19937 global_mt19937;

///Calling anything in this namespace is deprecated
namespace OptimTPrivate {    
inline uint32_t makeRandSeed() {
        static bool isFirstCalled=true;
        if(isFirstCalled) {
            uint32_t seed=std::time(nullptr),seed_=seed;
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
}   // OptimTPrivate


///uniform random number in range [0,1)
inline double randD() {
    static std::uniform_real_distribution<double> rnd(0,1);
    return rnd(global_mt19937);
}

///uniform random number in range [min,max]
inline double randD(const double min,const double max) {
    return (max-min)*randD()+min;
}

}   // OptimT

#endif // RANDOMS_HPP