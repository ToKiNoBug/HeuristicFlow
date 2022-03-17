// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DEF_TESTINGMATRIX_H
#define DEF_TESTINGMATRIX_H
#include <HeuristicFlow/SimpleMatrix>
void testAccess();

void testCopy();

void testCustomTypes();

void testLoop(unsigned int = 100);

void testInverse();

void testProduct();
#endif // DEF_TESTINGMATRIX_H
