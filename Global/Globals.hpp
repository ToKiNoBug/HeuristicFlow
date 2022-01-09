
#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <random>

namespace OptimT
{
///macro function for square
#define OT_square(x) (x)*(x)

///Empty class to put global variable and functions, instances of it is meanning less
class OtGlobal
{
public:
    ///global random device for OptimT
    static std::random_device global_random_device;
    ///global random device(mt19937) for OptimT
    static std::mt19937 global_mt19937;

    ///uniform random number in range [0,1)
    static double randD() {
        static std::uniform_real_distribution<double> rnd(0,1);
        return rnd(OtGlobal::global_mt19937);
        //return (std::rand()%RAND_MAX)/double(RAND_MAX);
    }
    ///uniform random number in range [min,max]
    static double randD(const double min,const double max) {
        return (max-min)*randD()+min;
    }
private:
    OtGlobal();
};


}

#endif // GLOBALS_HPP
