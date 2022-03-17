// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_PSOOPTION_HPP
#define Heu_PSOOPTION_HPP

#include "../../Global"

#include <stdint.h>

namespace Heu {

///Options to PSO
struct PSOOption
{
public:
    PSOOption() {
        populationSize=200;
        maxGeneration=300;
        maxFailTimes=100;
        inertiaFactor=0.8;
        learnFactorP=2;
        learnFactorG=2;
    }
    ///size of population
    size_t populationSize;
    ///maximum allowed generation
    size_t maxGeneration;
    ///maximun allowed failtimes;
    size_t maxFailTimes;
    ///speed factor
    double inertiaFactor;
    ///pBest factor
    double learnFactorP;
    ///gBest factor
    double learnFactorG;
};

}

#endif // PSOOPTION_HPP
