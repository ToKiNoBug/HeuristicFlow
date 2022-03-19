// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DEF_TESTINGGENETIC_H
#define DEF_TESTINGGENETIC_H

#include <Eigen/Dense>
#define Heu_USE_THREADS
#define Heu_NSGA_USE_THREADS
#include <HeuristicFlow/Genetic>


///Ackely function
void testAckley_withRecord();

///TSP problem
void testTSP_SOGA(const unsigned int pointNum=10);

void testDistance(const size_t=4);

///Binh and Korn function
void testNSGA2_Binh_and_Korn();

void testNSGA2_Kursawe();

void testNSGA2_ZDT3();

void testNSGA3_DTLZ7();

void searchPF();
#endif // DEF_TESTINGGENETIC_H
