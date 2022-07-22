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

/**
 * \brief The objective function
 *
 * f1(x,y)=4*(x^2+y^2);
 * f2(x,y)=(x-5)^2+(y-5)^2
 *
 * Non-liner constraints : (x-5)^2+y^2<=25 and 7.7<=(x-8)^2+(y+3)^2
 *
 * Here we use the penalty method to implement non-liner constraints.
 *
 * \tparam args_t Type of args
 * \param _x std array of 2 dim double vector [x,y]
 * \param f Fitness value to be computed
 */

void testNSGA2_Binh_and_Korn() {
  // 0<=x_0<=5,  0<=x_1<=3

  // this function requires 0<=x<=5 and 0<=y<=3. Use a non-square constraint to present it.
  using args_t = heu::ContinousBox<std::array<double, 2>, heu::BoxShape::RECTANGLE_BOX>;

  // the type of solver
  using solver_t = heu::NSGA2<
      std::array<double, 2>, 2, heu::FITNESS_LESS_BETTER, heu::RecordOption::DONT_RECORD_FITNESS,
      args_t,
      heu::GADefaults<std::array<double, 2>, args_t, heu::Std>::iFun,  // initializatoin functon
      nullptr,  // the fitness function can be assigned at runtime if nullptr is used
      heu::GADefaults<std::array<double, 2>, args_t, heu::Std>::cFunNd<>,  // otherwise the function
                                                                           // must be assigned in
                                                                           // the template
      heu::GADefaults<std::array<double, 2>, args_t, heu::Std>::mFun>;     // This is suitable for
                                                                           // iFun, fFun, cFun and
                                                                           // mFun.

  solver_t algo;

  // using Fitness_t = typename solver_t::Fitness_t;

  {  // Set the parameters
    heu::GAOption opt;
    opt.maxGenerations = 100;
    opt.populationSize = 200;
    // opt.maxFailTimes = 200;  // this member is useless to MOGA solvers

    opt.crossoverProb = 0.8;
    opt.mutateProb = 0.1;
    algo.setOption(opt);
  }

  {  // Set the box constraint
    args_t box;
    box.setRange({0, 0}, {5, 3});
    box.setDelta({0.05, 0.03});
    algo.setArgs(box);
  }

  // Set the fitness function at runtime.
  algo.setfFun(heu::testFunctions<std::array<double, 2>, Eigen::Array2d, args_t>::BinhKorn);

  // Initialize the population
  algo.initializePop();

  cout << "Start" << endl;
  std::clock_t t = std::clock();
  // Run the algorithm
  algo.run();
  t = std::clock() - t;
  cout << "Solving finished in " << double(t) / CLOCKS_PER_SEC << " seconds and "
       << algo.generation() << "generations" << endl;

  // output the pareto front
  cout << "paretoFront=[";
  for (const auto &i : algo.pfGenes()) {
    cout << i->_Fitness[0] << " , " << i->_Fitness[1] << ";\n";
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
  testNSGA2_Binh_and_Korn();

  // system("pause");
  return 0;
}
