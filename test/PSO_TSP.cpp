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
#include <iostream>
#include <ctime>
using namespace Eigen;
using namespace std;

#include <algorithm>

// this tests shows how to solve TSP problems with PSO
void testTSP_PSO(const int N) {
  static constexpr int SpaceDim = 2;
  // type of decision variable encoded in float
  using Var_t = Eigen::ArrayXd;
  //  type of args
  using DistanceMat_t = Eigen::ArrayXXd;
  //  type of solver
  using Solver_t = heu::PSO<Eigen::ArrayXd, heu::BoxShape::SQUARE_BOX, heu::FITNESS_LESS_BETTER,
                            heu::RECORD_FITNESS, DistanceMat_t>;

  using Args_t = Solver_t::Args_t;

  Solver_t solver;

  using sortUnit = std::pair<float, uint32_t>;
  //  Decode ArrayXd into a permulation through sorting.
  Solver_t::fFun_t fFun = [](const Var_t* x, const Args_t* args, double* fitness) {
    const size_t _N = (*args).rows();
    std::vector<sortUnit> sortSpace(_N);
    for (size_t i = 0; i < _N; i++) {
      sortSpace[i] = std::make_pair(x->operator[](i), i);
    }

    static const auto cmpFun = [](const sortUnit& a, const sortUnit& b) {
      return a.first < b.first;
    };

    std::sort(sortSpace.begin(), sortSpace.end(), cmpFun);

    *fitness = 0;

    for (size_t i = 0; i + 1 < _N; i++) {
      *fitness += (*args)(sortSpace[i].second, sortSpace[i + 1].second);
    }
  };

  DistanceMat_t dMat(N, N);

  dMat.setZero();

  {
    Eigen::Array<float, SpaceDim, Eigen::Dynamic> points;
    points.setRandom(SpaceDim, N);

    for (size_t r = 0; r < N; r++) {
      for (size_t c = 0; c < N; c++) {
        if (r == c) continue;
        dMat(r, c) = (points.col(r) - points.col(c)).square().sum();
      }
    }
  }

  // Set the option of PSO
  heu::PSOOption opt;
  opt.inertiaFactor = 0.8;
  opt.learnFactorG = 2;
  opt.learnFactorP = 2;
  opt.maxGeneration = 50 * N;
  opt.maxFailTimes = opt.maxGeneration / 10;
  opt.populationSize = 100;

  // Set the dimension of PSO
  solver.setDimensions(N);

  // Set posMin,posMax and velocityMax
  solver.setRange(0, 1);
  solver.setMaxVelocity(0.5);

  solver.setfFun(fFun);
  solver.setOption(opt);
  solver.setArgs(dMat);
  solver.initializePop();

  {
    double iniFitness;
    fFun(&solver.population().front().position, &solver.args(), &iniFitness);
    // The initial fitness value as reference
    cout << "iniFitness = " << iniFitness << endl;
  }

  clock_t c = clock();
  solver.run();
  c = clock() - c;

  cout << "finished in " << double(c) * 1000 / CLOCKS_PER_SEC << " miliseconds and "
       << solver.generation() << " generations" << endl;

  cout << "result fitness = " << solver.bestFitness() << endl;
  /*
  cout<<"Trainning Curve=[";
  for(auto i : solver.record()) {
      cout<<i<<" , ";
  }
  cout<<"];"<<endl;

  cout<<"Population condition:"<<endl;
  for(const auto & i : solver.population()) {
      cout<<"fitness="<<i.fitness<<" , pBest="<<i.pBest.fitness<<endl;
  }
  */
}

int main(int argc,char**argv) {
  bool is_auto=false;

  for(int i=0;i<argc;i++) {
    if(std::string_view{argv[i]}=="--auto") {
      is_auto=true;
      break;
    }
  }
  int NodeNum = 100;
  if(!is_auto) {
    cout << "Input node number : ";
    cin >> NodeNum;
  }
  testTSP_PSO(NodeNum);

  return 0;
}
