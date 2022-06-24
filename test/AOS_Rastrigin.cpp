#include <HeuristicFlow/AOS>
#include <HeuristicFlow/EAGlobal>

#include <iostream>

using std::cout, std::endl;

int main() {
  constexpr int Dim = 2;

  using Arg_t = int;
  heu::AOS<std::vector<double>, heu::BoxShape::RECTANGLE_BOX, heu::FITNESS_LESS_BETTER,
           heu::RECORD_FITNESS, Arg_t,
           heu::testFunctions<std::vector<double>, double, Arg_t>::rastrigin, false,
           heu::DivEncode<-5, 1>::code, heu::DivEncode<5, 1>::code, heu::DivEncode<3, 2>::code>
      solver;

  heu::AOSOption opt;
  opt.electronNum = 50;
  opt.maxEarlyStop = 15;
  opt.maxGeneration = 50;
  opt.photonRate = 0.1;
  opt.maxLayerNum = 5;

  solver.setOption(opt);
  /*
solver.setDimensions(Dim);
solver.setMax(5);
solver.setMin(-5);
solver.setLearnRate(1.0);
*/
  solver.setMin({-5, -5});
  solver.setMax({5, 5});
  solver.setLearnRate({1, 1});

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