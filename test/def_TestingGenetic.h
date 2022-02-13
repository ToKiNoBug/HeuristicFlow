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

#ifndef DEF_TESTINGGENETIC_H
#define DEF_TESTINGGENETIC_H

#include <Eigen/Dense>
#define Heu_USE_THREADS
#define Heu_NSGA_USE_THREADS
#include <HeuristicFlow/Genetic>


///Using libEigen with HeuristicFlow
void testWithEigenLib();

///Ackely function
void testAckley_withRecord();

///TSP problem
void testTSP(const unsigned int pointNum=10);

void testDistance(const size_t=4);

///Binh and Korn function
void testNSGA2_Binh_and_Korn();

void testNSGA2_Kursawe();

void testNSGA2_ZDT3();

void testNSGA3_DTLZ7();
#endif // DEF_TESTINGGENETIC_H
