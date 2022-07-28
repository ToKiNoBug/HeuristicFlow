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

#include <HeuristicFlow/AOS>
#include <HeuristicFlow/EAGlobal>

#include <iostream>

using std::cout, std::endl;

int main() {
  constexpr int Dim = 3;

  heu::AOS<heu::FixedContinousBox17<Eigen::Array<double, Dim, 1>, heu::encode(-5), heu::encode(5),
                                    heu::encode(1.5)>,
           heu::FITNESS_LESS_BETTER, heu::RECORD_FITNESS, void,
           heu::testFunctions<Eigen::Array<double, Dim, 1>>::rastrigin>
      solver;

  // solver.setRange(-5, 5);
  // solver.setDelta(1.5);

  heu::AOSOption opt;
  opt.electronNum = 50;
  opt.maxEarlyStop = 20;
  opt.maxGeneration = 50;
  opt.photonRate = 0.1;
  opt.maxLayerNum = 5;

  solver.setOption(opt);

  solver.initializePop();

  /*
cout << "box range = " << solver.min() << " , " << solver.max() << " , " << solver.learnRate()
     << endl;

cout << "InitialS=[";
for (const auto& X : solver.electrons()) {
  cout << X.state.transpose() << ";\n";
}
cout << "]';\n\n\n\n" << endl;
*/
  std::clock_t clocks = std::clock();
  solver.run();
  clocks = std::clock() - clocks;

  /*cout << "AOS finished with " << double(clocks) / CLOCKS_PER_SEC * 1e6 / solver.generation()
       << " ms per generation" << endl;*/

  // cout << "result Var_t = [" << solver.bestElectron().state << "];\n\n";

  cout << "result fitness = " << solver.bestElectron().energy << " in " << solver.generation()
       << " generations" << endl;

  /*
cout << "record=[";
for (auto i : solver.record()) {
  cout << i << ',';
}
cout << "];" << endl;
*/
}