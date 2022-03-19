// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DEF_TESTINGPSO_H
#define DEF_TESTINGPSO_H

#include <Eigen/Dense>
#include <iostream>
#define Heu_USE_THREADS
#include <HeuristicFlow/PSO>

void testPSOBase();

void testRastriginFun();

void testTSP_PSO(const size_t N);

#endif // DEF_TESTINGPSO_H
