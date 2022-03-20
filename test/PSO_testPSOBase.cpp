// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Eigen/Dense>
#include <HeuristicFlow/PSO>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;

void testPSOBase() {
    using Var_t = std::vector<double>;

    internal::PSOAbstract<Var_t,
            double,
            DONT_RECORD_FITNESS,void,nullptr,nullptr> * abstractNoRec=nullptr;

    internal::PSOAbstract<Var_t,
            double,
            RecordOption::RECORD_FITNESS,void,nullptr,nullptr> * abstractDoRec=nullptr;

    internal::PSOBase<Var_t,10,
            double,
            RecordOption::DONT_RECORD_FITNESS,void,nullptr,nullptr> * baseNoRec=nullptr;

    internal::PSOBase<Var_t,0,
            double,
            RecordOption::RECORD_FITNESS,void,nullptr,nullptr> * baseDoRec=nullptr;

    //DoRec is derived from NoRec
    abstractNoRec=abstractDoRec;

    //base is derived from abstract
    abstractNoRec=baseNoRec;
    abstractDoRec=baseDoRec;

}


int main()
{
    testPSOBase();
    system("pause");
    return 0;
}