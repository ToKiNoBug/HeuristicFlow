
#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <random>
#include <ctime>
#include "./LogisticChaos.h"
#ifdef OptimT_DO_PARALLELIZE
#include <omp.h>
#endif

namespace OptimT
{
///macro function for square
#define OT_square(x) (x)*(x)

///global random device for OptimT
//extern std::random_device global_random_device;
///global random device(mt19937) for OptimT
extern std::mt19937 global_mt19937;

extern LogisticChaos global_logistic;
///Infinet value for float
const float pinfF=1.0f/0.0f;

///Infinet value for double
const double pinfD=1.0/0.0;

///negative infinet value for floa
const float nInfF=-pinfF;

///negative infinet value for double
const double ninfD=-pinfD;

///Empty class to put global variable and functions, instances of it is meanning less
class OtGlobal
{
public:
    /// [Deprecated]
    ///uniform random number in range [0,1)
    static double randD() {
        static std::uniform_real_distribution<double> rnd(0,1);
        return rnd(global_mt19937);
        //return (std::rand()%RAND_MAX)/double(RAND_MAX);
    }
    /// [Deprecated]
    ///uniform random number in range [min,max]
    static double randD(const double min,const double max) {
        return (max-min)*randD()+min;
    }
    static uint32_t makeRandSeed() {
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
private:
    OtGlobal();
};

///uniform random number in range [0,1)
inline double randD() {
    static std::uniform_real_distribution<double> rnd(0,1);
    return rnd(global_mt19937);
}

///uniform random number in range [min,max]
inline double randD(const double min,const double max) {
    return (max-min)*randD()+min;
}

///logistic random number in range (0,1)
inline double logisticD() {
    return global_logistic();
}

#define OPTIMT_MAKE_GLOBAL \
std::mt19937 OptimT::global_mt19937(OptimT::OtGlobal::makeRandSeed()); \
OptimT::LogisticChaos OptimT::global_logistic(randD());

}

#endif // GLOBALS_HPP
