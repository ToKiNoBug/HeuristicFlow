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

void testNSGA3_DTLZ7() {
  // Dimension of decison variable
  static constexpr size_t N = 30;
  // Number of objectives. DTLZ series are greatly extentable, you can simply edit the value of M
  // and N to see how it behaves for higher dimensions.
  static constexpr size_t M = 3;

  // Note that NSGA3 supports only FITNESS_LESS_BETTER
  using solver_t = heu::NSGA3<Eigen::Array<double, N, 1>, M,
                              // FITNESS_LESS_BETTER,
                              heu::DONT_RECORD_FITNESS, heu::SINGLE_LAYER, void>;

  using Var_t = Eigen::Array<double, N, 1>;
  // using Fitness_t = solver_t::Fitness_t;

  auto iFun = heu::GADefaults<Var_t>::iFunNd<>;

  auto cFun = heu::GADefaults<Var_t>::cFunNd<heu::DivEncode<1, 10>::code>;

  auto mFun = [](const Var_t* src, Var_t* v) {
    *v = *src;
    const size_t idx = heu::randIdx(v->size());
    double& p = v->operator[](idx);
    p += 0.4 * heu::randD(-1, 1);
    p = std::min(p, 1.0);
    p = std::max(p, 0.0);
  };

  heu::GAOption opt;
  opt.maxGenerations = 1000;
  opt.populationSize = 400;
  cout << "maxGenerations=";
  cin >> opt.maxGenerations;

  // opt.maxFailTimes = -1; // this member is useless to MOGA solvers
  cout << "populationSize=";
  cin >> opt.populationSize;
  opt.crossoverProb = 0.8;
  opt.mutateProb = 0.1;

  solver_t solver;
  // solver.setObjectiveNum(M);
  solver.setiFun(iFun);
  // solver.setfFun(heu::testFunctions<Var_t, Eigen::Array<double, M, 1>>::template DTLZ4<>);
  solver.setfFun(heu::testFunctions<Var_t, Eigen::Array<double, M, 1>>::DTLZ7);
  //  solver.setfFun(DTLZ1<M, N>);
  solver.setcFun(cFun);
  solver.setmFun(mFun);
  solver.setOption(opt);
  solver.setReferencePointPrecision(12);
  // solver.setReferencePointPrecision(5, 7);
  cout << "RPCount=" << solver.referencePointCount() << endl;
  solver.initializePop();
  cout << "Population initialized." << endl;
  // cout << "RP=[" << solver.referencePoints() << "]';\n\n\n" << endl;

  clock_t c = clock();
  solver.run();
  c = clock() - c;

  cout << "solving finished in " << c << " ms with " << solver.generation() << " generations."
       << endl;

  cout << "PFV=[";
  for (const auto& i : solver.pfGenes()) {
    cout << i->_Fitness.transpose() << ";\n";
  }
  cout << "];\n\n\n" << endl;
}

// This function goes through the whole high-dimensional searching space to find an accurate
// solution of PF. It's really slow but you can find if the solution of NSGA3 is correct.
template <size_t M, size_t N, size_t precision>
void searchPFfun(size_t nIdx, Eigen::Array<double, N, 1>* var,
                 vector<pair<Eigen::Array<double, M, 1>, size_t>>* dst) {
  static_assert(precision >= 2, "You should assign at least 2 on a single dim");
  if (nIdx + 1 >= N) {
    Eigen::Array<double, M, 1> f;
    for (size_t i = 0; i <= precision; i++) {
      var->operator[](nIdx) = double(i) / precision;
      heu::testFunctions<Eigen::Array<double, N, 1>, Eigen::Array<double, M, 1>>::DTLZ7(var, &f);
      dst->emplace_back(make_pair(f, 0ULL));
    }
    return;
  }

  for (size_t i = 0; i <= precision; i++) {
    var->operator[](nIdx) = double(i) / precision;
    searchPFfun<M, N, precision>(nIdx + 1, var, dst);
  }
}

int main() {
  testNSGA3_DTLZ7();
  system("pause");
  return 0;
}
