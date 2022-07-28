#include <HeuristicFlow/Global>

using namespace std;
using namespace heu;

void printBin(const uint64_t bin, bool noSpace = false) {
  char str[67] = "";

  uint64_t mask = 1ULL << 63;

  for (int idx = 0; idx < 64; idx++) {
    str[idx] = (mask & bin) == 0 ? '0' : '1';
    mask = mask >> 1;
  }

  for (int idx = 0; idx < 64; idx++) {
    printf("%c", str[idx]);
    if (!noSpace) {
      if (idx == 0 || idx == 11) {
        printf("%c", ' ');
      }
    }
  }

  printf("\n");
}

inline void printBin(const double fl, bool noSpace = false) {
  printBin(reinterpret_cast<const uint64_t &>(fl), noSpace);
}

int main() {
  printf("%s\n", "Correct val of 1.52e-320 : ");
  printBin(1.52e-320);

  printf("%s\n", "computed val of 1.52e-320 : ");
  printBin(encode(1.52e-320));

  printf("%s\n", "Correct val of 3.14e-100 : ");
  printBin(3.14e-100);

  printf("%s\n", "computed code of of 3.14e-100 : ");
  printBin(encode(3.14e-100));

  printf("%s\n", "Correct val of 1.507 : ");
  printBin(1.507);

  printf("%s\n", "computed code of of 1.507 : ");
  printBin(encode(1.507));

  printf("%s\n", "Correct val of 3.14 : ");
  printBin(3.14);

  printf("%s\n", "computed code of of 3.14 : ");
  printBin(encode(3.14));

  printf("%s\n", "Correct val of 3.14e100 : ");
  printBin(3.14e100);

  printf("%s\n", "computed code of of 3.14e100 : ");
  printBin(encode(3.14e100));
}