#include <HeuristicFlow/SimpleMatrix>

#include <iostream>

#include <string>

void printBin(size_t v);

void printBinInInterval(const size_t* v, const size_t bits, const size_t bitPerSpace);

int main() {
  using namespace heu;
  using std::cout, std::endl;

  multiBitSet<7> vec(25);

  // vec.resize(25);

  // cout << "(0b00100ULL << char(-1)) = " << (0b00100ULL << char(1)) << endl;

  for (int byteId = 0; byteId < vec.blocks() * sizeof(decltype(vec)::block_t); byteId++) {
    reinterpret_cast<uint8_t*>(vec.data())[byteId] = heu::randIdx(0, 256);
  }

  cout << "size = " << vec.size() << " , capacity = " << vec.capacity() << endl;

  cout << "storage space in binary is : \n";
  for (int blockId = 0; blockId < vec.blocks(); blockId++) printBin(vec.data()[blockId]);

  cout << "\n\n\nBinary in elements is : \n";

  printBinInInterval(vec.data(), 7 * 25, 7);

  cout << "\n\n\n";

  cout << "The elements are : [";

  for (auto val : vec) {
    cout << int(val) << ", ";
  }

  cout << "\b\b];" << endl;

  /*
cout << "\n[9] = " << int(vec[9]) << endl;

const auto val = 0b1111111;

cout << "set [9] to " << val << endl;

vec[9] = val;

*/

  /*
#warning set value to 1 doesn't works well
  cout << "\n\nstorage space in binary is : \n";
  for (int blockId = 0; blockId < vec.blocks(); blockId++) printBin(vec.data()[blockId]);

  cout << "\n\n\n";

  cout << "The elements are : [";

  for (auto val : vec) {
    cout << int(val) << ", ";
  }

  cout << "\b\b];" << endl;

  */

  return 0;
}

void printBin(size_t v) {
  size_t mask = size_t(1) << (8 * sizeof(size_t) - 1);

  char buf[8 * sizeof(size_t) + 1] = "";

  for (int bitIdx = 0; bitIdx < sizeof(size_t) * 8; bitIdx++) {
    buf[bitIdx] = ((mask & v) > 0) ? ('1') : ('0');
    mask = mask >> 1;
  }

  printf("%s ", buf);
}

void printBinInInterval(const size_t* v, const size_t bits, const size_t bitPerSpace) {
  const size_t size_tNum = std::ceil(float(bits) / sizeof(size_t) / 8.0);

  std::string str;
  str.reserve(bits * 2);
  size_t bitIdx = 0;
  for (int size_tIdx = 0; size_tIdx < size_tNum; size_tIdx++) {
    size_t mask = 1ULL << 63;
    const size_t val = v[size_tIdx];

    for (int bit = 0; bit < 64; bit++) {
      str.push_back((mask & val) ? '1' : '0');
      bitIdx++;
      mask = mask >> 1;

      if (bitIdx % bitPerSpace == 0) {
        str.push_back(' ');
      }
    }
  }

  std::cout << str << std::endl;
}