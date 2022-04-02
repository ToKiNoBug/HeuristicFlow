/*
 Copyright © 2021-2022  TokiNoBug
This file is part of HeuristicFlow.

    HeuristicFlow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HeuristicFlow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HeuristicFlow.  If not, see <https://www.gnu.org/licenses/>.

*/
#include <Eigen/Dense>
#include <HeuristicFlow/EAGlobal>

#include <iostream>
using namespace Eigen;
using namespace std;

// this example shows how to use box-constraint types
void test_Box() {
  //[0,1]^50, mutate step = 0.02
  heu::BoxNdS<50, heu::ContainerOption::Eigen, true, heu::DivEncode<0, 1>::code, heu::DivEncode<1, 1>::code,
              heu::DivEncode<1, 50>::code>
      box0;

  box0.min();

  cout << box0.dimensions() << endl;
  cout << box0.learnRate() << endl;
  cout << "sizeof(box0) = " << sizeof(box0)
       << endl;  // every thing about this boxes are known at compile time, so its size is 1

  // Dynamic dim non-suqare box
  heu::BoxXdN<heu::ContainerOption::Eigen> box;

  // Set the min and max by its non-const-ref to max and min members.
  box.max().setConstant(50, 1, 1.0);
  box.min().setConstant(50, 1, 0.0);

  // 10 dimensional binary box
  heu::BoxNb<10> BNb;
  // Dynamic dimensional binary box
  heu::BoxXb<> BXb;

  BXb.setDimensions(400);

  BNb.max();
  BNb.min();

  BXb.max();
  BXb.min();

  cout << BNb.dimensions();
  cout << BXb.dimensions();

  cout << "sizeof BNb=" << sizeof(BNb) << endl;  // Every information of BNb is fixed at compile time so its size is 1
  cout << "sizeof BXb=" << sizeof(BXb) << endl;  // BXb 's dimensions is known at runtime so its size is 8 (size_t)
}

int main() {
  test_Box();
  system("pause");
  return 0;
}
