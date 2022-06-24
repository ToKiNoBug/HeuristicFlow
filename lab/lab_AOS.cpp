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
constexpr int electronNum = 50;
constexpr int maxGeneration = 50;
constexpr int maxLayerNum = 5;

constexpr double photonRate = 0.5;

constexpr double varMin = -5;
constexpr double varMax = 5;

/*
 * Referenced paper : Mahdi Azizi.Atomic orbital search: A novel metaheuristic algorithm[J].Applied
 * Mathematical Modelling.2021,93:657-693. Link : https://doi.org/10.1016/j.apm.2020.12.021 *
 *
 */

inline double gaussianCurve(const double x, const double mu = 0.0,
                            const double sigma = 1.0) noexcept {
  return 1 / (sigma * std::sqrt(2 * M_PI)) *
         std::exp(-heu::square(x - mu) / (2 * heu::square(sigma)));
}

void ackely(const Eigen::Array2d* _x, double* f) noexcept {
  heu::testFunctions<Eigen::Array2d>::rastrigin(_x, f);
  return;

  double x = _x->operator[](0), y = _x->operator[](1);
  *f = -20 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) -
       exp(0.5 * (cos(M_2_PI * x) + cos(M_2_PI * y))) + 20 + M_E;
  //*f = std::log10(*f);
}

class Electron {
 public:
  Electron() { isComputed = false; }
  Eigen::Array2d state;
  double energy;
  bool isComputed;
  void setUncomputed() noexcept { isComputed = false; }
};

template <typename pointer_t>
bool electronSortCompareFun(pointer_t A, pointer_t B) {
  return (A->energy < B->energy);
}

class Layer : public std::vector<Electron*> {
 public:
  Eigen::Array2d bindingState;
  double bindingEnergy;
  int LEidx;

  const Eigen::Array2d& LEState() const { return this->at(LEidx)->state; }
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
  std::list<Electron> electrons;
  std::vector<Layer> atom;
  Eigen::Array2d bindingState;
  double bindingEnergy;
  std::list<Electron>::iterator lowestEnergyIterator;
  int generations;
  std::vector<double> record;

  void initializePop() {
    electrons.resize(electronNum);
    for (auto& elec : electrons) {
      elec.setUncomputed();
      heu::randD(elec.state.data(), elec.state.size(), varMin, varMax);
    }
    atom.resize(electronNum);
    for (auto& layer : atom) {
      layer.reserve(electronNum);
      layer.clear();
    }
    atom.clear();
    record.clear();
    record.reserve(maxGeneration + 1);
  }

  void run() {
    generations = 0;
    while (true) {
      computeAll();

      updateAtomBSBE();

      record.emplace_back(this->lowestEnergyIterator->energy);

      makeLayers();

      if (generations >= maxGeneration) {
        break;
      }

      updateLayeBSBE();

      updateElectrons();

      generations++;
    }
  }

 protected:
  void computeAll() {
    for (auto& e : electrons) {
      if (e.isComputed) {
        continue;
      }
      ackely(&e.state, &e.energy);
      e.isComputed = true;
    }
  }

  void updateAtomBSBE() {
    bindingState = electrons.front().state;
    bindingEnergy = electrons.front().energy;
    lowestEnergyIterator = electrons.begin();

    for (auto it = electrons.begin(); it != electrons.end(); ++it) {
      bindingState += it->state;
      bindingEnergy += it->energy;
      if (it->energy < lowestEnergyIterator->energy) {
        lowestEnergyIterator = it;
      }
    }

    bindingEnergy /= electrons.size();
    bindingState /= electrons.size();
  }
  void makeLayers() {
    std::vector<std::list<Electron>::iterator> elecSortSpace;
    elecSortSpace.reserve(electrons.size());
    for (auto it = electrons.begin(); it != electrons.end(); ++it) {
      elecSortSpace.emplace_back(it);
    }

    std::sort(elecSortSpace.begin(), elecSortSpace.end(),
              electronSortCompareFun<std::list<Electron>::iterator>);

    while (elecSortSpace.size() > electronNum) {
      electrons.erase(elecSortSpace.back());
      elecSortSpace.pop_back();
    }

    //////////////////////////////////////////////////
    atom.clear();
    int layerNum = heu::randIdx(1, maxLayerNum + 1);

    Eigen::ArrayXd numOfEachLayer(layerNum);

    numOfEachLayer.setZero();

    for (int idx = 0; idx < layerNum; idx++) {
      numOfEachLayer[idx] = gaussianCurve(idx + 1, 0, layerNum / 6);
    }

    numOfEachLayer /= numOfEachLayer.sum();
    numOfEachLayer *= electrons.size();
    Eigen::ArrayXi intNumOfEachLayer = numOfEachLayer.round().cast<int>();

    for (int idx = 1; idx < intNumOfEachLayer.size(); idx++) {
      intNumOfEachLayer[idx] += intNumOfEachLayer[idx - 1];
    }

    for (auto& val : intNumOfEachLayer) {
      if (val > electrons.size()) val = electrons.size();
    }

    atom.resize(layerNum);
    for (auto& i : atom) {
      i.clear();
    }
    int curLayerIdx = 0;

    for (int countedElectrons = 0; countedElectrons < electrons.size(); countedElectrons++) {
      if (countedElectrons >= intNumOfEachLayer[curLayerIdx]) {
        curLayerIdx++;
      }

      atom[curLayerIdx].emplace_back(&*elecSortSpace[countedElectrons]);
    }

    while (atom.back().size() <= 0) {
      atom.pop_back();
    }
  }

  inline void updateLayeBSBE() {
    for (auto& l : atom) {
      l.updateLayerBSBE();
    }
  }

  void updateElectrons() {
    for (int layerIdx = 0; layerIdx < atom.size(); layerIdx++) {
      for (auto elecPtr : atom[layerIdx]) {
        electrons.emplace_back();
        Electron* newElecPtr = &electrons.back();
        if (heu::randD() > photonRate) {
          applyPhoton(atom[layerIdx], layerIdx + 1, *elecPtr, newElecPtr);
        } else {
          applyNonPhoton(*elecPtr, newElecPtr);
        }

        for (double& var : newElecPtr->state) {
          var = std::max(var, varMin);
          var = std::min(var, varMax);
        }

        newElecPtr->setUncomputed();
      }
    }
  }

 private:
  void applyPhoton(const Layer& layer, const int layerK, const Electron& prevElec,
                   Electron* newElec) {
    Eigen::Array2d alpha, beta, gamma;

    heu::randD(alpha.data(), alpha.size());
    heu::randD(beta.data(), beta.size());
    heu::randD(gamma.data(), gamma.size());

    if (prevElec.energy >= layer.bindingEnergy) {
      auto deltaX = alpha * (beta * this->lowestEnergyIterator->state - gamma * this->bindingState);
      newElec->state = prevElec.state + deltaX / double(layerK);
    } else {
      auto deltaX = alpha * (beta * layer.LEState() - gamma * layer.bindingState);
      newElec->state = prevElec.state + deltaX;
    }
    newElec->setUncomputed();
  }

  void applyNonPhoton(const Electron& prevElec, Electron* newElec) {
    Eigen::Array2d rand;
    heu::randD(rand.data(), rand.size(), 0, 1);
    newElec->state = prevElec.state + rand;
    newElec->setUncomputed();
  }
};

int main() {
  // test tst;

  labAOS solver;
  solver.initializePop();
  solver.run();
  /*
cout << "result=[" << solver.lowestEnergyIterator->state << "];\n\n";

cout << "Trainning Curve:\nfitness=[";
for (auto fitness : solver.record) {
  cout << fitness << ',';
}
cout << "];\n" << endl;
*/

  cout << "best fitness=" << solver.lowestEnergyIterator->energy << endl;

  // system("pause");
}