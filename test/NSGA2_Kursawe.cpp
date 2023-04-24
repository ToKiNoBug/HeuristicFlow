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
#include <HeuristicFlow/Genetic>
#include <iostream>
#include <ctime>
using namespace Eigen;
using namespace std;

// This test shows how to use NSGA2 solvers without args.

/*
Kursawe function is also a testing function for multi-objective evoluntionary algorithms
The Kursawe function has 1~3 dim(s) input and 2 dim(s) output.

*/

void testNSGA2_Kursawe() {
  // 1<=i<=3,  -5<=x_i<=5

  // All functions at compile time is their default value (nullptr), so they must be assigned at
  // runtime
  heu::NSGA2<std::array<double, 3>, 2, heu::FITNESS_LESS_BETTER, heu::DONT_RECORD_FITNESS> algo;

  // initialize function
  auto iFun = [](std::array<double, 3> *x) {
    for (auto &i : *x) {
      i = heu::randD(-5, 5);
    }
  };

  // The kursawe function is the fitness function
  auto fFun = heu::testFunctions<std::array<double, 3>, Eigen::Array<double, 2, 1>>::Kursawe;

  // crossover function
  auto cFun = [](const std::array<double, 3> *p1, const std::array<double, 3> *p2,
                 std::array<double, 3> *ch1, std::array<double, 3> *ch2) {
    for (int i = 0; i < 3; i++) {
      static const double r = 0.2;
      ch1->operator[](i) = r * p1->operator[](i) + (1 - r) * p2->operator[](i);
      ch2->operator[](i) = r * p2->operator[](i) + (1 - r) * p1->operator[](i);
    }
  };

  // mutate function
  auto mFun = [](const std::array<double, 3> *src, std::array<double, 3> *x) {
    *x = *src;
    const size_t idx = heu::randIdx(3);
    x->operator[](idx) += 0.1 * heu::randD(-1, 1);
    x->operator[](idx) = std::min(x->operator[](idx), 5.0);
    x->operator[](idx) = std::max(x->operator[](idx), -5.0);
  };

  heu::GAOption opt;
  opt.maxGenerations = 2000;
  opt.populationSize = 600;
  // opt.maxFailTimes = opt.maxGenerations;  // this member is useless to MOGA solvers

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
  cout << "Solving finished in " << double(t) / CLOCKS_PER_SEC << " seconds and "
       << algo.generation() << " generations" << endl;
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
  return 0;
}
