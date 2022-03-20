// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef Heu_GLOBALS_HPP
#define Heu_GLOBALS_HPP

#include <random>
#include <ctime>
#include <thread>
#include <type_traits>
#include <vector>
#include "Constants.hpp"
#include "Randoms.hpp"
#include "HeuMaths.hpp"

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


#endif // GLOBALS_HPP
