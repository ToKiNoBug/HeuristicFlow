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

#include <iostream>
#include <array>
#include <quadmath.h>

using std::endl, std::cout;

int main() {
  __float128 A1 = 0.5 * sinq(1.0) - 2 * cosq(1.0) + 1 * sinq(2.0) - 1.5 * cosq(2.0);
  __float128 A2 = 1.5 * sinq(1.0) - 1 * cosq(1.0) + 2 * sinq(2.0) - 0.5 * cosq(2.0);

  std::array<char, 1024> buf;

  quadmath_snprintf(buf.data(), buf.size(), "%3.128Qf", A1);
  cout << "A1 = " << buf.data() << endl;
  quadmath_snprintf(buf.data(), buf.size(), "%3.128Qf", A2);
  cout << "A2 = " << buf.data() << endl;

  return 0;
}