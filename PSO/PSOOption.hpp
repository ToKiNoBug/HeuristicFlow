/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef PSOOPTION_HPP
#define PSOOPTION_HPP

#include "../OptimTemplates/Global"

#include <stdint.h>

namespace OptimT {

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

    size_t populationSize;
    size_t maxGeneration;
    size_t maxFailTimes;
    double inertiaFactor;
    double learnFactorP;
    double learnFactorG;
};

}

#endif // PSOOPTION_HPP
