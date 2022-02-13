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

#ifndef DEF_TESTINGPSO_H
#define DEF_TESTINGPSO_H

#include <Eigen/Dense>
#include <iostream>
#define Heu_USE_THREADS
#include <HeuristicFlow/PSO>

void testPSOBase();

void testRastriginFun();

void testTSP(const size_t N);

#endif // DEF_TESTINGPSO_H
