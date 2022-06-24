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