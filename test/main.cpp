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
#include "def_TestingFuns.h"

using namespace std;

///initialize std::rand()
void initializeSrand();

int main()
{
    initializeSrand();
    testTSP(20);
    return 0;
}

void initializeSrand() {
    std::time_t t=std::time(nullptr);
    std::srand((t>>32)^(t&0xFFFFFFFF));
}
