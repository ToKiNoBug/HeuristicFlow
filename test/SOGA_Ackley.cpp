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

// Ackely function is a single objective testing function with great number of local minimum points,
// but its global minimum point is [0,0] and corresponding value is 0. A good solver should be able
// to find the global minimum point.
void testAckley_withRecord(bool is_auto) {
  using args_t = heu::ContinousBox<array<double, 2>, heu::BoxShape::SQUARE_BOX>;
  constexpr heu::SelectMethod sm = heu::SelectMethod::RunTimeSelectMethod;

  using solver_t = heu::SOGA<array<double, 2>, heu::FITNESS_LESS_BETTER, heu::RECORD_FITNESS, sm,
                             args_t, heu::GADefaults<array<double, 2>, args_t>::iFun<>, nullptr,
                             heu::GADefaults<array<double, 2>, args_t>::cFunNd,
                             heu::GADefaults<array<double, 2>, args_t>::mFun<>>;
  solver_t algo;

  heu::GAOption opt;
  opt.populationSize = 50;
  opt.maxFailTimes =
      50;  // The solver will end if it hasn't been finding a better solution for 50 generations.
  opt.maxGenerations = 100;

  {
    uint32_t sm = heu::SelectMethod::EliteReserved;
    if(!is_auto) {
      cout << "Please choose a selection method :\n";
      for (heu::SelectMethod i = heu::SelectMethod(0); i < heu::SelectMethod::RunTimeSelectMethod;
           i = heu::SelectMethod(int(i) + 1)) {
        cout << i << " for " << heu::Enum2String(i) << ",\n";
      }

      cout << "Please nput the number of selection method : ";

      std::cin >> sm;
    }

    algo.setSelectMethod(heu::SelectMethod(sm));
  }

  // if constexpr (sm == heu::SelectMethod::Tournament)
  algo.setTournamentSize(3);

  // if constexpr (sm == heu::SelectMethod::LinearRank)
  algo.setLinearSelectProbability(0.2, 0.8);

  // if constexpr (sm == heu::SelectMethod::ExponentialRank)
  algo.setExponetialSelectBase(0.8);

  // if constexpr (sm == heu::SelectMethod::Boltzmann)
  algo.setBoltzmannSelectStrength(-10);

  // if constexpr (sm == heu::SelectMethod::EliteReserved)
  algo.setEliteNum(3);

  algo.setOption(opt);

  {
    args_t args;
    args.setRange(-5, 5);
    args.setDelta(0.05);
    algo.setArgs(args);
  }

  algo.setfFun(
      // Ackely function
      heu::testFunctions<array<double, 2>, double, args_t>::ackley
      /*
      [](const array<double, 2>* _x, const solver_t::ArgsType*, double* f) {
        double x = _x->operator[](0), y = _x->operator[](1);
        *f = -20 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) - exp(0.5 * (cos(M_2_PI * x) + cos(M_2_PI
      * y))) + 20 + M_E;
      }*/

  );

  algo.initializePop();

  std::clock_t t = std::clock();
  algo.run();
  t = std::clock() - t;
  // cout<<algo.bestFitness();

  cout << "Solving spend " << algo.generation() << " generations in " << double(t) / CLOCKS_PER_SEC
       << " sec\n";
  cout << "Result = [" << algo.result()[0] << " , " << algo.result()[1] << "]\n";

  cout << "Fitness history :\n";

  for (auto i : algo.record()) {
    cout << i << '\n';
  }
  cout << endl;
}

int main(int argc,char**argv) {
  bool is_auto=false;

  for(int i=0;i<argc;i++) {
    if(std::string_view{argv[i]}=="--auto") {
      is_auto=true;
      break;
    }
  }

  testAckley_withRecord(is_auto);
  // system("pause");
  return 0;
}
