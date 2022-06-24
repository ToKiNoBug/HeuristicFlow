#include <HeuristicFlow/AOS>
#include <HeuristicFlow/EAGlobal>

#include <iostream>

using std::cout, std::endl;

int main() {
  constexpr int Dim = 2;

  heu::AOS<Eigen::Array<double, Dim, 1>, heu::BoxShape::SQUARE_BOX, heu::FITNESS_LESS_BETTER,
           heu::RECORD_FITNESS, void, heu::testFunctions<Eigen::Array<double, Dim, 1>>::rastrigin,
           true, heu::DivEncode<-5, 1>::code, heu::DivEncode<5, 1>::code,
           heu::DivEncode<3, 2>::code>
      solver;

  heu::AOSOption opt;
  opt.electronNum = 50;
  opt.maxEarlyStop = -1;
  opt.maxGeneration = 1000;
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

  cout << "AOS finished with " << double(clocks) / CLOCKS_PER_SEC * 1e6 / solver.generation()
       << " ms per generation" << endl;

  // cout << "result Var_t = [" << solver.bestElectron().state << "];\n\n";
  /*
cout << "result fitness = " << solver.bestElectron().energy << " in " << solver.generation()
     << " generations" << endl;
     */
  /*
cout << "record=[";
for (auto i : solver.record()) {
  cout << i << ',';
}
cout << "];" << endl;
*/
}