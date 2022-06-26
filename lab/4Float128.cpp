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