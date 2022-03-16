/*
 Copyright © 2022  TokiNoBug
This file is part of Heuristic.

    Heuristic is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Heuristic is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Heuristic.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef DEF_TESTINGBOXES_H
#define DEF_TESTINGBOXES_H


#include <Eigen/Dense>
#define Heu_USE_THREADS
#define Heu_NSGA_USE_THREADS
#include <HeuristicFlow/EAGlobal>

void test_Box_double();

void test_Box_bool();

#endif // DEF_TESTINGBOXES_H
