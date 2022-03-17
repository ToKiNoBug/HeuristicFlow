// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <ctime>
#include "def_TestingGenetic.h"
#include "def_TestingPSO.h"
#include "def_TestingMatrix.h"
#include "def_TestingBoxes.h"

using namespace std;

int main() {
    //test_Box_double();
    testAckley_withRecord();
    //testNSGA2_Binh_and_Korn();
    //testNSGA3_DTLZ7();
    //searchPF();

    system("pause");
    return 0;
}
