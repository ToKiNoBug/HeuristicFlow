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

#ifndef DEF_TESTINGGENETIC_H
#define DEF_TESTINGGENETIC_H

void testSingleNumber();

///Using libEigen with OptimTemplates
void testWithEigenLib();

///Ackely function
void testAckley_withRecord();

///TSP problem
void testTSP(const unsigned int pointNum=10);
#endif // DEF_TESTINGGENETIC_H