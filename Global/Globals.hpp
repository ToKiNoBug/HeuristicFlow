
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

#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <random>
#include <ctime>
#include <thread>
#include "./Enumerations.hpp"
#include "./Chaotic.hpp"
#include "./Randoms.hpp"
#include "./OptimTMaths.hpp"
#include <type_traits>
#ifdef OptimT_DO_PARALLELIZE
#include <omp.h>
#include <thread>
#endif

namespace OptimT
{
///macro function for square
#define OT_square(x) (x)*(x)

///Size identifier for dynamic size (fitness or var)
const size_t Dynamic = 0;

///Infinet value for float
const float pinfF=1.0f/0.0f;

///Infinet value for double
const double pinfD=1.0/0.0;

///negative infinet value for floa
const float nInfF=-pinfF;

///negative infinet value for double
const double ninfD=-pinfD;

///Empty class to put global variable and functions, instances of it is meanningless
class OtGlobal
{
public:
    static uint32_t threadNum() {
        return concurrency;
    }
private:
    OtGlobal();
    static uint32_t concurrency;
};


#define OptimT_MAKE_GLOBAL \
std::mt19937 OptimT::global_mt19937(OptimT::OptimTPrivate::makeRandSeed()); \
OptimT::LogisticChaos OptimT::global_logistic(randD()); \
uint32_t OptimT::OtGlobal::concurrency=std::thread::hardware_concurrency();

}

#endif // GLOBALS_HPP
