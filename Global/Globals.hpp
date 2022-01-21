
#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <random>
#include <ctime>
#include <thread>
#include "./LogisticChaos.hpp"

#ifdef OptimT_DO_PARALLELIZE
#include <omp.h>
#include <thread>
#endif

namespace OptimT
{
///macro function for square
#define OT_square(x) (x)*(x)

///whether to record trainning curve of not
enum RecordOption : uint8_t {
    RECORD_FITNESS=true,
    DONT_RECORD_FITNESS=false
};
///optimization direction
enum FitnessOption : uint8_t {
    FITNESS_LESS_BETTER=false,
    FITNESS_GREATER_BETTER=true,
};

enum ConstraintOption : uint8_t {
    NONCONSTRAINT,
    IS_CONSTRAINT
};

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
    static uint32_t threadNum() {
        return concurrency;
    }

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
    static uint32_t concurrency;
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

#define OptimT_MAKE_GLOBAL \
std::mt19937 OptimT::global_mt19937(OptimT::OtGlobal::makeRandSeed()); \
OptimT::LogisticChaos OptimT::global_logistic(randD()); \
uint32_t OptimT::OtGlobal::concurrency=std::thread::hardware_concurrency();

}

#endif // GLOBALS_HPP
