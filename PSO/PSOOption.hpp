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

#ifndef OptimT_PSOOPTION_HPP
#define OptimT_PSOOPTION_HPP

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
