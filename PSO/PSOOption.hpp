#ifndef PSOOPTION_HPP
#define PSOOPTION_HPP

#include "../OptimTemplates/Global"

#include <stdint.h>

namespace OptimT {

struct PSOOption
{
public:
    PSOOption() {
        populationSize=100;
        maxGeneration=300;
        maxFailTimes=100;
        inertiaFactor=1;
        learnFactorP=2;
        learnFactorG=2;
    }

    size_t populationSize;
    size_t maxGeneration;
    size_t maxFailTimes;
    double inertiaFactor;
    double learnFactorP;
    double learnFactorG;
};

}

#endif // PSOOPTION_HPP
