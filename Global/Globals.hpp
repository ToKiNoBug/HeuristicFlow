
#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <random>

namespace OptimT
{
///macro function for square
#define OT_square(x) (x)*(x)



///global random device for OptimT
extern std::random_device global_random_device;
///global random device(mt19937) for OptimT
extern std::mt19937 global_mt19937;

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
    ///uniform random number in range [0,1)
    static double randD() {
        static std::uniform_real_distribution<double> rnd(0,1);
        return rnd(global_mt19937);
        //return (std::rand()%RAND_MAX)/double(RAND_MAX);
    }
    ///uniform random number in range [min,max]
    static double randD(const double min,const double max) {
        return (max-min)*randD()+min;
    }
private:
    OtGlobal();
};


#define OPTIMT_MAKE_GLOBAL \
std::random_device OptimT::global_random_device; \
std::mt19937 OptimT::global_mt19937(OptimT::global_random_device());


}

#endif // GLOBALS_HPP
