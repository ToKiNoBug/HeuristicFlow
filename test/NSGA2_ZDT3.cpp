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

/// Zitzler–Deb–Thiele's function N. 3 is a MO testing function
void testNSGA2_ZDT3() {
  // 0<=x_i<=1, 1<=i<=30
  static const size_t XNum = 30;
  // static const double r = 0.2;

  heu::NSGA2<std::array<double, XNum>, 2, heu::FITNESS_LESS_BETTER, heu::RECORD_FITNESS> algo;

  void (*iFun)(std::array<double, XNum>*) = [](std::array<double, XNum>* x) {
    for (size_t i = 0; i < XNum; i++) {
      x->operator[](i) = heu::randD();
    }
  };

  void (*fFun)(const std::array<double, XNum>* x, Eigen::Array<double, 2, 1>* f) =
      [](const std::array<double, XNum>* x, Eigen::Array<double, 2, 1>* f) {
        f->operator[](0) = x->operator[](0);
        const double&& f1 = std::move(f->operator[](0));
        double g = 0;
        for (size_t i = 1; i < XNum; i++) {
          g += x->operator[](i);
        }
        g = 1 + g * 9.0 / (XNum - 1);
        double&& f1_div_g = f1 / g;
        double&& h = 1 - std::sqrt(f1_div_g) - f1_div_g * std::sin(10 * M_PI * f1);
        f->operator[](1) = g * h;
      };

  auto cFun = heu::GADefaults<std::array<double, XNum>>::cFunNd<(heu::DivEncode<1, 2>::code)>;

  void (*mFun)(const std::array<double, XNum>* src, std::array<double, XNum>*) =
      [](const std::array<double, XNum>* src, std::array<double, XNum>* x) {
        *x = *src;
        const size_t mutateIdx = heu::randIdx(XNum);

        x->operator[](mutateIdx) += 0.005 * heu::randD(-1, 1);

        x->operator[](mutateIdx) = std::min(x->operator[](mutateIdx), 1.0);
        x->operator[](mutateIdx) = std::max(x->operator[](mutateIdx), 0.0);
      };

  heu::GAOption opt;
  opt.populationSize = 200;
  // opt.maxFailTimes = 500;  // this member is useless to MOGA solvers
  opt.maxGenerations = 2000;

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
  cout << "Solving finished in " << double(t) / CLOCKS_PER_SEC << "seconds and "
       << algo.generation() << "generations" << endl;
  std::vector<Eigen::Array<double, 2, 1>> paretoFront;
  algo.paretoFront(paretoFront);
  cout << "paretoFront=[";
  for (const auto& i : paretoFront) {
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
  testNSGA2_ZDT3();
  system("pause");
  return 0;
}
