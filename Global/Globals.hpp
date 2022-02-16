
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

#ifndef Heu_GLOBALS_HPP
#define Heu_GLOBALS_HPP

#include <random>
#include <ctime>
#include <thread>
#include <type_traits>
#include <vector>
#include "./Constants.hpp"
#include "./Chaotic.hpp"
#include "./Randoms.hpp"
#include "./HeuMaths.hpp"

#include <array>

#ifndef EIGEN_CORE_H    //Detects whether libEigen is included
#ifdef Heu_GENETIC_USE_EIGEN     //If user hopes to use Eigen without including it, report an error
#error You must include Eigen before you define Heu_GENETIC_USE_EIGEN! Include Eigen before Heu.
#endif
#endif

#ifdef Heu_USE_THREADS
#include <omp.h>
#include <thread>
#endif

namespace Heu
{
///Empty class to put global variable and functions, instances of it is meanningless
class HfGlobal
{
public:
    inline static uint32_t threadNum() {
        return concurrency;
    }
private:
    HfGlobal() {};
    ~HfGlobal() {};
    static uint32_t concurrency;
};


#define Heu_MAKE_GLOBAL \
std::mt19937 Heu::global_mt19937(Heu::HeuPrivate::makeRandSeed()); \
Heu::LogisticChaos Heu::global_logistic(randD()); \
uint32_t Heu::HfGlobal::concurrency=std::thread::hardware_concurrency();

}

#endif // GLOBALS_HPP
