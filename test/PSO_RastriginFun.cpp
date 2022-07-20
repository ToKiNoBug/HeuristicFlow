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
#include <HeuristicFlow/PSO>
#include <HeuristicFlow/EAGlobal>
#include <iostream>
#include <ctime>
using namespace Eigen;
// using namespace std;
using namespace heu;

using std::cout, std::endl;

// Test NSGA2 with rastrigin function
void testRastriginFun() {
  static constexpr size_t N = 20;
  using Var_t = Eigen::Array<double, N, 1>;

  using solver_t = heu::PSO<Var_t, BoxShape::SQUARE_BOX, heu::FITNESS_LESS_BETTER,
                            heu::RECORD_FITNESS, void, heu::testFunctions<Var_t>::rastrigin>;

  heu::PSOOption opt;
  opt.populationSize = 400;
  opt.maxGeneration = 50 * N;
  opt.maxFailTimes = -1;
  // opt.maxGeneration/10;
  opt.inertiaFactor = 0.8;
  opt.learnFactorG = 2;
  opt.learnFactorP = 2;

  solver_t solver;

  solver.setRange(-5.12, 5.12);
  solver.posMax() = 5.12;  //  you can also set the min and max pos in this way
  solver.setMaxVelocity(0.1);

  solver.setOption(opt);

  solver.initializePop();

  // solver.initialize(iFun,fFun,opt);

  clock_t time = clock();
  solver.run();
  time = clock() - time;

  cout << "finished in " << double(time) * 1000 / CLOCKS_PER_SEC << " miliseconds and "
       << solver.generation() << " generations" << endl;

  cout << "result fitness = " << solver.bestFitness() << endl;
  /*
  const Var_t &result = solver.globalBest().position;
  cout<<"result = [";
  for(auto i : result) {
      cout<<i<<" , ";
  }
  cout<<"];"<<endl;


  cout<<"Trainning Curve=[";
  for(auto i : solver.record()) {
      cout<<i<<" , ";
  }
  cout<<"];"<<endl;

  */

  /*
  cout<<"Population condition:"<<endl;
  for(const auto & i : solver.population()) {
      cout<<"fitness="<<i.fitness<<" , pBest="
         <<i.pBest.fitness<<" , position="<<i.position.transpose()<<" ,
  velocity="<<i.velocity.transpose()<<endl;
  }
  */
}

int main() {
  testRastriginFun();
  system("pause");
  return 0;
}
