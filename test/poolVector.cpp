#include <HeuristicFlow/SimpleMatrix>

#include <iostream>
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

int main() {
  cout << __LINE__ << endl;
  heu::poolVector<testScalar> vec;
  cout << __LINE__ << endl;
}