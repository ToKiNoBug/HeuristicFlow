// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Eigen/Dense>
#include <HeuristicFlow/Genetic>
#include <iostream>
#include <ctime>
using namespace Eigen;
using namespace std;

void testNSGA2_Kursawe() {
  // 1<=i<=3,  -5<=x_i<=5
  NSGA2<std::array<double, 3>, 2, FITNESS_LESS_BETTER, DONT_RECORD_FITNESS> algo;
  auto iFun = [](std::array<double, 3> *x) {
    for (auto &i : *x) {
      i = ei_randD(-5, 5);
    }
  };
  auto fFun = [](const std::array<double, 3> *x, Eigen::Array<double, 2, 1> *f) {
    double f1 = 0, f2 = 0;
    for (int i = 0; i < 2; i++) {
      f1 += -10 *
            exp(-0.2 * sqrt((x->operator[](i)) * (x->operator[](i)) + (x->operator[](i + 1)) * (x->operator[](i + 1))));
    }
    for (int i = 0; i < 3; i++) {
      f2 += pow(abs(x->operator[](i)), 0.8) + 5 * sin(x->operator[](i) * x->operator[](i) * x->operator[](i));
    }
    f->operator[](0) = f1;
    f->operator[](1) = f2;
  };

  auto cFun = [](const std::array<double, 3> *p1, const std::array<double, 3> *p2, std::array<double, 3> *ch1,
                 std::array<double, 3> *ch2) {
    for (int i = 0; i < 3; i++) {
      static const double r = 0.2;
      ch1->operator[](i) = r * p1->operator[](i) + (1 - r) * p2->operator[](i);
      ch2->operator[](i) = r * p2->operator[](i) + (1 - r) * p1->operator[](i);
    }
  };

  auto mFun = [](const std::array<double, 3> *src, std::array<double, 3> *x) {
    *x = *src;
    const size_t idx = ei_randIdx(3);
    x->operator[](idx) += 0.1 * ei_randD(-1, 1);
    x->operator[](idx) = std::min(x->operator[](idx), 5.0);
    x->operator[](idx) = std::max(x->operator[](idx), -5.0);
  };

  GAOption opt;
  opt.maxGenerations = 2000;
  opt.populationSize = 600;
  opt.maxFailTimes = -1;

  algo.setiFun(iFun);
  algo.setmFun(mFun);
  algo.setfFun(fFun);
  algo.setcFun(cFun);
  algo.setOption(opt);
  algo.initializePop();

  cout << "Start" << endl;
  std::clock_t t = std::clock();
  algo.run();
  t = std::clock() - t;
  cout << "Solving finished in " << double(t) / CLOCKS_PER_SEC << " seconds and " << algo.generation() << " generations"
       << endl;
  std::vector<Eigen::Array<double, 2, 1>> paretoFront;
  algo.paretoFront(paretoFront);
  cout << "paretoFront=[";
  for (const auto &i : paretoFront) {
    cout << i[0] << " , " << i[1] << ";\n";
  }
  cout << "];" << endl;
  /*
  cout<<"\n\n\n population=[";
  for(const auto & i : algo.population()) {
      cout<<i.fitness()[0]<<" , "<<i.fitness()[1]<<";\n";
  }
  cout<<"];"<<endl;
  */
}

int main() {
  testNSGA2_Kursawe();
  system("pause");
  return 0;
}