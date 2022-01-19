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
#include "def_TestingGenetic.h"

#include <unordered_set>

using namespace std;

OptimT_MAKE_GLOBAL

void tryOMP();

int main()
{
    testNSGA2_Binh_and_Korn();
    system("pause");
    return 0;
}
#include <string>
#include <iomanip>
void tryOMP() {
    const char str[]="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
#pragma omp parallel for
    for(uint32_t i=0;i<sizeof(str)-1;i++) {
        cout<<char(str[i]);
    }

}
