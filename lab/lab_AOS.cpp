#include <Eigen/Dense>

#include <HeuristicFlow/Global>
#include <HeuristicFlow/EAGlobal>
#include <HeuristicFlow/SimpleMatrix>

#include <list>
#include <iostream>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;

constexpr int dimensions = 2;
constexpr int electronNum = 100;

inline double gaussianCurve(const double x, const double mu = 0.0, const double sigma = 1.0) {
  return 1 / (sigma * std::sqrt(2 * M_PI)) * std::exp(-heu::square(x - mu) / (2 * heu::square(sigma)));
}

void ackely(const Eigen::Array2d* _x, double* f) {
  double x = _x->operator[](0), y = _x->operator[](1);
  *f = -20 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) - exp(0.5 * (cos(M_2_PI * x) + cos(M_2_PI * y))) + 20 + M_E;
}

struct Electron {
  Eigen::Array2d state;
  double energy;
};

class test {
 public:
  test() { cout << "Constructor" << endl; }
  ~test() { cout << "Destructor" << endl; }
};

bool electronSortCompareFun(const Electron* A, const Electron* B) { return (A->energy < B->energy); }

class Layer : public std::vector<Electron*> {
 public:
  Eigen::Array2d bindingState;
  double bindingEnergy;
  int LEidx;

  void updateLayerBSBE() {
    LEidx = 0;
    bindingState = this->front()->state;
    bindingEnergy = this->front()->energy;
    for (int idx = 1; idx < this->size(); idx++) {
      bindingState += this->operator[](idx)->state;
      bindingEnergy += this->operator[](idx)->energy;
      if (this->operator[](idx)->energy < this->operator[](LEidx)->energy) {
        LEidx = idx;
      }
    }
    bindingState /= this->size();
    bindingEnergy /= this->size();
  };
};

class labAOS {
 public:
  std::vector<Electron> electrons;
  std::vector<Layer> atom;
  Eigen::Array2d bindingState;
  double bindingEnergy;
  int LEidx;

  void updateAtomBSBE() {
    bindingState = electrons.front().state;
    bindingEnergy = electrons.front().energy;
    LEidx = 0;
    for (int idx = 1; idx < electrons.size(); idx++) {
      bindingState += electrons[idx].state;
      bindingEnergy += electrons[idx].energy;
      if (electrons[idx].energy < electrons[LEidx].energy) {
        LEidx = idx;
      }
    }

    bindingEnergy /= electrons.size();
    bindingState /= electrons.size();
  }

  void updateLayeBSBE() {
    for (auto& l : atom) {
      l.updateLayerBSBE();
    }
  }

  void makeLayers() {
    atom.reserve(electronNum);
    atom.clear();
    int layerNum = heu::randIdx(3, int(electrons.size()));
    Eigen::ArrayXd numOfEachLayer(layerNum);
    numOfEachLayer.setZero();
    for (int idx = 0; idx < layerNum; idx++) {
      numOfEachLayer[idx] = gaussianCurve(idx + 1, 0, layerNum / 6);
    }

    numOfEachLayer /= numOfEachLayer.sum();
    numOfEachLayer *= electronNum;
    Eigen::ArrayXi intNumOfEachLayer = numOfEachLayer.round().cast<int>();

    for (int idx = 1; idx < intNumOfEachLayer.size(); idx++) {
      intNumOfEachLayer[idx] += intNumOfEachLayer[idx - 1];
    }

    for (auto& val : numOfEachLayer) {
      if (val > electronNum) val = electronNum;
    }

    std::vector<Electron*> elecSortSpace(electronNum);
    for (int idx = 0; idx < electronNum; idx++) {
      elecSortSpace[idx] = this->electrons.data() + idx;
    }

    std::sort(elecSortSpace.begin(), elecSortSpace.end(), electronSortCompareFun);

    atom.resize(layerNum);
    for (auto& i : atom) {
      i.clear();
    }
    int curLayerIdx = 0;

    for (int countedElectrons = 0; countedElectrons < electronNum; countedElectrons++) {
      if (countedElectrons >= intNumOfEachLayer[curLayerIdx]) {
        curLayerIdx++;
      }

      atom[curLayerIdx].emplace_back(elecSortSpace[countedElectrons]);
    }

    while (atom.back().size() <= 0) {
      atom.pop_back();
    }
  }

  void initializePop() { electrons.resize(electronNum); }
};

int main() {
  //
  std::vector<test> testArr;
  cout << __LINE__ << endl;
  testArr.reserve(10);
  cout << __LINE__ << endl;
  testArr.resize(10);
  cout << __LINE__ << endl;
  testArr.resize(0);
  cout << __LINE__ << endl;

  // test tst;

  // system("pause");
}