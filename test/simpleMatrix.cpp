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