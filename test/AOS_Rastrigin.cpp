#include <HeuristicFlow/AOS>
#include <HeuristicFlow/EAGlobal>

#include <iostream>

using std::cout, std::endl;

int main() {
  constexpr int Dim = 2;

  heu::AOS<std::vector<double>, heu::BoxShape::SQUARE_BOX, heu::FITNESS_LESS_BETTER,
           heu::RECORD_FITNESS, void, heu::testFunctions<std::vector<double>>::rastrigin, true,
           heu::DivEncode<-5, 1>::code, heu::DivEncode<5, 1>::code, heu::DivEncode<3, 2>::code>
      solver;

  heu::AOSOption opt;
  opt.electronNum = 50;
  opt.maxEarlyStop = 15;
  opt.maxGeneration = 50;
  opt.photonRate = 0.1;
  opt.maxLayerNum = 5;

  solver.setOption(opt);

  solver.setDimensions(Dim);

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

  solver.run();

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