#include <HeuristicFlow/Global>
#include <iostream>

using namespace heu;

int main() {
  constexpr int minV = min_v<int, 1, 2, 3, 4, -114514, 6, 7, 8, 9, 5>;

  std::cout << "minV = " << minV << std::endl;
  constexpr int maxV = max_v<int, 1, 2, 1919810, 4, 5, 6, 7, 8, 9, 3>;
  std::cout << "maxV = " << maxV << std::endl;

  return 0;
}