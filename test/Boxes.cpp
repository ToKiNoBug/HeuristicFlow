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
#include "Eigen/src/Core/Array.h"
using namespace Eigen;
using namespace std;

using heu::BoxShape;

// this example shows how to use box-constraint types
void test_Box() {
  cout << "Eigen::Matrix4d is EigenClass ? " << heu::array_traits<Eigen::Matrix4d>::isEigenClass
       << endl;
  cout << "Eigen::Array44d is EigenClass ? " << heu::array_traits<Eigen::Array44d>::isEigenClass
       << endl;

  cout << "std::vector is EigenClass ? " << heu::array_traits<std::vector<int>>::isEigenClass
       << endl;
  //[0,1]^50, mutate step = 0.02

  cout << "Testing Box1 : box4dS" << endl;

  heu::ContinousBox<Eigen::Array4d, heu::BoxShape::SQUARE_BOX> box0;

  box0.setRange(0, 1);
  box0.setDelta(0.6);
  box0.min();

  cout << box0.dimensions() << endl;
  // cout << box0.learnRate() << endl;
  cout << "sizeof(box0) = " << sizeof(box0)
       << endl;  // every thing about this boxes are known at compile time, so its size is 1
  Eigen::Array4d v;
  v.fill(0.5);
  cout << "v before constraint : \n" << v.transpose() << endl;
  box0.applyConstraint(&v);

  cout << "v after constraint : \n" << v.transpose() << endl;

  box0.applyDelta(&v);

  cout << "v after delta : \n" << v.transpose() << endl;

  cout << "\n\n\n\n\nTesting Box1 : boxXiS" << endl;

  // Dynamic dim square integer box
  heu::DiscretBox<std::vector<int>, BoxShape::SQUARE_BOX> box1;

  box1.setDimensions(20);

  box1.setRange(-3, 4);

  std::vector<int> a;

  box1.initialize(&a);

  cout << "initialized : =[";
  for (int i : a) {
    cout << i << ',';
  }
  cout << "\b]" << endl;

  cout << "\n\n\n\n\nTesting Box2 : boxXXfNS" << endl;
  // Dynamic non-square mat flaot box
  heu::ContinousBox<Eigen::ArrayXXf, heu::BoxShape::RECTANGLE_BOX> box2;

  box2.setDimensions(3, 6);

  box2.min().fill(-2);

  box2.max().fill(6);

  for (int idx = 0; idx < box2.dimensions(); idx++) {
    box2.min()(idx) -= (idx % 5);
    box2.max()(idx) += (idx % 5);
  }

  box2.delta() = (box2.max() - box2.min()) / 10;

  cout << "box2 shape : [" << box2.boxRows() << " , " << box2.boxCols() << "]" << endl;

  cout << "box2.min() = \n" << box2.min() << endl;
  cout << "box2.max() = \n" << box2.max() << endl;
  cout << "box2.delta() = \n" << box2.delta() << endl;

  decltype(box2)::Var_t mat2;

  box2.initialize(&mat2);

  cout << "initialized : \n" << mat2 << endl;

  box2.applyConstraint(&mat2);
  box2.applyDelta(&mat2);

  cout << " mat after delta : \n" << mat2 << endl;

  cout << "\n\n\n\n\nTesting Box3 : box10b" << endl;
  // 10 dimensional binary box

  heu::DiscretBox<std::array<bool, 10>> box3;

  cout << "sizeof box10b = " << sizeof(box3) << endl;

  std::array<bool, 10> b;

  box3.initialize(&b);

  cout << "Intialized : [";
  for (bool i : b) {
    cout << i << ',';
  }
  cout << "\b];" << endl;

  box3.applyDelta(&b);

  cout << "After delta : [";
  for (bool i : b) {
    cout << i << ',';
  }
  cout << "\b];" << endl;

  cout << "\n\n\n\n\nTesting Box4 : inf boxXXf" << endl;
  heu::GaussianBox<Eigen::ArrayXXd> box4;

  box4.setDimensions(3, 4);
  box4.setMu(10);
  box4.setSigma(1e8);
  box4.setDelta(1e7);

  Eigen::ArrayXXd mat4;

  box4.initialize(&mat4);

  // mat4.matrix().determinant();

  cout << "initialized : \n" << mat4 << endl;

  box4.applyDelta(&mat4);

  cout << "after delta : \n" << mat4 << endl;

  cout << "\n\n\n\n\nTesting Box5 : Fixed box23i (in range [-3,6])" << endl;
  heu::FixedDiscreteBox<Eigen::Array<int, 2, 3>, -3, -2> box5;

  cout << "size of box5 = " << sizeof(box5) << endl;

  cout << "box5.dimensions() = " << box5.dimensions() << endl;

  cout << "box5.min() = " << box5.min() << endl;
  cout << "box5.max() = " << box5.max() << endl;

  Eigen::Array<int, 2, 3> mat5;

  box5.initialize(&mat5);

  cout << "initialized by box5 : \n" << mat5 << endl;

  cout << "\n\n\n\n\nTesting box6 : Fixed BoxXXf (in range [-10,5], delta = 0.5)" << endl;

  heu::FixedContinousBox17<Eigen::ArrayXXd, heu::encode(-10), heu::encode(5), heu::encode(0.5)>
      box6;
  box6.setDimensions(3, 4);

  Eigen::ArrayXXd mat6;

  cout << "box6.min() = " << box6.min() << endl;
  cout << "box6.max() = " << box6.max() << endl;
  cout << "box6.delta() = " << box6.delta() << endl;

  cout << "size of box6 = " << sizeof(box6) << endl;

  box6.initialize(&mat6);

  cout << "initialized : \n" << mat6 << endl;

  box6.applyDelta(&mat6);

  cout << "after delta : \n" << mat6 << endl;
}

int main() {
  test_Box();
  // system("pause");
  return 0;
}
