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

#ifndef Heu_GAOPTION_HPP
#define Heu_GAOPTION_HPP
#include <stdint.h>
namespace Heu {

///options about GA algorithm
struct GAOption
{
public:
    GAOption() {
        populationSize=100;
        maxFailTimes=50;
        maxGenerations=300;
        crossoverProb=0.8;
        mutateProb=0.05;
    }
    size_t populationSize;
    size_t maxFailTimes;
    size_t maxGenerations;
    double crossoverProb;
    double mutateProb;
};

}

#endif // GAOPTION_HPP
