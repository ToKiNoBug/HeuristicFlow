// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "def_TestingBoxes.h"
#include <iostream>
using namespace Heu;
using namespace std;

void test_Box_double() {
    BoxNdS<50,DoubleVectorOption::Eigen,
            10,encode<0,1>::code,encode<1,1>::code,encode<1,50>::code> box0;
    //50 dim square box in [0,1]
    BoxXdN<DoubleVectorOption::Eigen> box;

    //cout<<box.Flag<<endl;

    box.max().setConstant(50,1,1.0);
    box.min().setConstant(50,1,0.0);

    box0.min();

    cout<<box0.dimensions()<<endl;
    cout<<box0.learnRate()<<endl;
    cout<<"sizeof(box0) = "<<sizeof(box0)<<endl;

}

void test_Box_bool() {
    BoxNb<10> BNb;
    BoxXb<> BXb;

    BXb.setDimensions(400);

    BNb.max();
    BNb.min();

    BXb.max();
    BXb.min();

    cout<<BNb.dimensions();
    cout<<BXb.dimensions();

    cout<<"sizeof BNb="<<sizeof(BNb)<<endl;
    cout<<"sizeof BXb="<<sizeof(BXb)<<endl;

}


void test_random() 
{
    size_t N=20;
    for(size_t i=0;i<N;i++) {
        cout<<randD(-1,1)<<endl;
    }
}

void test_inf()
{
    Eigen::Array2d vec;
    vec[0]=pinfD;
    vec[1]=pinfD;

    cout<<vec.sum()<<endl;
}