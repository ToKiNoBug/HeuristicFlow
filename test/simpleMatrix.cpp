#include <HeuristicFlow/SimpleMatrix>
#include <string>
#include <iostream>
using namespace heu;
using namespace std;

int main() {
  MatrixDynamicSize<std::string> matXX(3, 4);

  MatrixFixedSize<std::string, 3, 4> mat34;

  mat34.fill("rua!");

  matXX.fill("rue~");
  cout << "Before Copying : " << endl;
  cout << "mat34(0,0) = " << mat34(0, 0) << endl;
  cout << "matXX(0,0) = " << matXX(0, 0) << endl;

  // matXX = mat34;
  mat34 = matXX;

  cout << "After Copying : " << endl;
  cout << "mat34(0,0) = " << mat34(0, 0) << endl;
  cout << "matXX(0,0) = " << matXX(0, 0) << endl;

  return 0;
}