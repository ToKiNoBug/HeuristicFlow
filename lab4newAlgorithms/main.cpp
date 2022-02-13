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

#include <iostream>
#include "includes.h"
#include "lab4NSGA2.h"
#include "lab4NSGA3.h"

Heu_MAKE_GLOBAL
using namespace std;
using namespace Eigen;
using namespace Heu;

int main() {
    static const double threshold=10000000000000000ULL;
    static const uint64_t res=Heu::amplifier<12345>::result;
    static const double decode=res/threshold;


    static const bool exceedsTh=(res>=threshold);


    static const bool excedds10Th=((res)>=10ULL*threshold);

    static const PowCode pc=powEncode<-12345678,255>::code;

    static const double val=powDecode<pc>::real;


    return 0;

    testNSGA3Expri();


    return 0;
}
