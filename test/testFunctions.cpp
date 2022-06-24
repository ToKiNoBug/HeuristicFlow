#include <HeuristicFlow/EAGlobal>

#include <Eigen/Dense>

#include <iostream>
using std::cout, std::endl;

int main() {
  constexpr double xMin = -512;
  constexpr double xMax = 512;
  constexpr double yMin = -512;
  constexpr double yMax = 512;

  constexpr int XYNum = 250;

  Eigen::ArrayXd x, y;
  Eigen::ArrayXXd z;

  x.setLinSpaced(XYNum, xMin, xMax);
  y = x;

  z.setZero(x.size(), y.size());

  for (int r = 0; r < x.size(); r++) {
    for (int c = 0; c < y.size(); c++) {
      Eigen::Array2d temp(x(r), y(c));
      heu::testFunctions<Eigen::Array2d>::eggHolder(&temp, &z(r, c));
    }
  }

  cout << "x=linspace(" << xMin << " , " << xMax << " , " << XYNum << ");" << endl;
  cout << "y=linspace(" << yMin << " , " << yMax << " , " << XYNum << ");" << endl;

  cout << "z=[";
  for (int r = 0; r < x.size(); r++) {
    cout << "[";
    for (int c = 0; c < y.size(); c++) {
      cout << z(r, c) << ',';
    }
    cout << "\b];";
  }
  cout << "\b];\n\n\n\n" << endl;

  cout << "surf(x,y,z,'EdgeColor','none','FaceAlpha',0.9)" << endl;

  // system("pause");
}