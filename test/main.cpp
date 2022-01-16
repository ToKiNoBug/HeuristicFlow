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

#include <iostream>
#include <ctime>
#include "testNsga2.h"

#include "def_TestingGenetic.h"

#include <unordered_set>

using namespace std;

OPTIMT_MAKE_GLOBAL

int main()
{
    auto t=OptimT::TentChaosI<8>::max();
    std::cout<<t<<std::endl;
    system("pause");

    testNSGA2_ZDT3();
    system("pause");
    return 0;
}
