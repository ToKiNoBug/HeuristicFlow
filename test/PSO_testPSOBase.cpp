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
#include <Eigen/Dense>
#include <HeuristicFlow/PSO>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;

// this tests shows the hierarchy of PSO solvers.
void testPSOBase() {
  using Var_t = std::vector<double>;

  heu::internal::PSOAbstract<Var_t, double, heu::DONT_RECORD_FITNESS, void, nullptr, nullptr>* abstractNoRec = nullptr;

  heu::internal::PSOAbstract<Var_t, double, heu::RecordOption::RECORD_FITNESS, void, nullptr, nullptr>* abstractDoRec =
      nullptr;

  heu::internal::PSOBase<Var_t, 10, double, heu::RecordOption::DONT_RECORD_FITNESS, void, nullptr, nullptr>* baseNoRec =
      nullptr;

  heu::internal::PSOBase<Var_t, 0, double, heu::RecordOption::RECORD_FITNESS, void, nullptr, nullptr>* baseDoRec =
      nullptr;

  // base is derived from abstract
  abstractNoRec = baseNoRec;
  abstractDoRec = baseDoRec;

  // DoRec is derived from NoRec
  abstractNoRec = abstractDoRec;
}

int main() {
  testPSOBase();
  system("pause");
  return 0;
}
