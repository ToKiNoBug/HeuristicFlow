#include <HeuristicFlow/SimpleMatrix>
#include <HeuristicFlow/Global>
#include <cmath>
#include <iostream>
#include <vector>
using std::cout;
using std::endl;

class testScalar {
 public:
  testScalar() { cout << "Constructor" << endl; }
  testScalar(const testScalar &) { cout << "copyConstructor by const reference" << endl; }
  testScalar(testScalar &&) { cout << "copyConstructor by right value reference" << endl; }
  ~testScalar() { cout << "Destructor" << endl; }
  double value;
};

class counterScalar {
 public:
  counterScalar() { createdNum++; }
  counterScalar(const counterScalar &) { createdNum++; }
  // counterScalar(counterScalar &&) { createdNum++; }
  ~counterScalar() { createdNum--; }
  static int createdNum;
  char val[3];
};

int counterScalar::createdNum = 0;

int main() {
  constexpr int loopN = 1e4;
  const int mallocMax = std::pow(2, 12);
  auto vec = new heu::poolVector<counterScalar>;
  // auto vec = new std::vector<counterScalar>;
  //   heu::poolVector<counterScalar> vec;
  cout << "createdNum=" << counterScalar::createdNum << endl;
  for (int loop = 0; loop < loopN; loop++) {
    // vec->push_back();
    vec->reserve(heu::randIdx(mallocMax));
    //  vec->clear();
  }
  delete vec;
  cout << "createdNum=" << counterScalar::createdNum << endl;
}