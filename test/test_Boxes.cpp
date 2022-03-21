// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Eigen/Dense>
#include <HeuristicFlow/EAGlobal>

#include <iostream>
using namespace Eigen;
using namespace std;

void test_Box() {
  BoxNdS<50, DoubleVectorOption::Eigen, true, DivEncode<0, 1>::code, DivEncode<1, 1>::code, DivEncode<1, 50>::code>
      box0;
  // 50 dim square box in [0,1]
  BoxXdN<DoubleVectorOption::Eigen> box;

  // cout<<box.Flag<<endl;

  box.max().setConstant(50, 1, 1.0);
  box.min().setConstant(50, 1, 0.0);

  box0.min();

  cout << box0.dimensions() << endl;
  cout << box0.learnRate() << endl;
  cout << "sizeof(box0) = " << sizeof(box0) << endl;

  BoxNb<10> BNb;
  BoxXb<> BXb;

  BXb.setDimensions(400);

  BNb.max();
  BNb.min();

  BXb.max();
  BXb.min();

  cout << BNb.dimensions();
  cout << BXb.dimensions();

  cout << "sizeof BNb=" << sizeof(BNb) << endl;
  cout << "sizeof BXb=" << sizeof(BXb) << endl;
}

int main() {
  test_Box();
  system("pause");
  return 0;
}