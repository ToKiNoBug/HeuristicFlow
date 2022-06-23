#include <HeuristicFlow/AOS>
#include <HeuristicFlow/EAGlobal>

#include <iostream>

using std::cout, std::endl;

int main() {
  heu::AOS<Eigen::Array4d, heu::BoxShape::SQUARE_BOX, heu::FITNESS_LESS_BETTER, heu::RECORD_FITNESS,
           void, heu::testFunctions<Eigen::Array4d>::rastrigin, true, heu::DivEncode<-5, 1>::code,
           heu::DivEncode<5, 1>::code, heu::DivEncode<3, 2>::code>
      solver;

  heu::AOSOption opt;
  opt.electronNum = 50;
  opt.maxEarlyStop = -1;
  opt.maxGeneration = 100;
  opt.photonRate = 0.5;

  solver.setOption(opt);

  solver.initializePop();

  solver.run();

  cout << "result Var_t = [" << solver.bestElectron().state << "];\n\n";

  cout << "result fitness = " << solver.bestElectron().energy << endl;

  cout << "record=[";
  for (auto i : solver.record()) {
    cout << i << ',';
  }
  cout << "];" << endl;
}