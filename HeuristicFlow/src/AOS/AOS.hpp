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

#ifndef HEU_AOS_HPP
#define HEU_AOS_HPP

#include <algorithm>

#include <HeuristicFlow/EAGlobal>

#include "InternalHeaderCheck.hpp"
#include "AOSBase.hpp"
#include "AOSDefaultElectron.hpp"

namespace heu {
template <typename Var_t, BoxShape BS, FitnessOption fOpt = FitnessOption::FITNESS_LESS_BETTER,
          RecordOption rOpt = RecordOption::DONT_RECORD_FITNESS, class Arg_t = void,
          typename internal::AOSParameterPack<Var_t, double, Arg_t>::fFun_t _fFun_ = nullptr,
          bool isFixedRange = false, DivCode MinCT = DivEncode<0, 1>::code,
          DivCode MaxCT = DivEncode<1, 1>::code, DivCode LRCT = DivEncode<1, 5>::code,
          typename internal::AOSParameterPack<Var_t, double, Arg_t>::iFun_t _iFun_ =
              internal::AOSParameterPack<Var_t, double,
                                         Arg_t>::defaultInitializeFunctionThatShouldNotBeCalled>
class AOS : public internal::AOSBase<
                Var_t, double, Arg_t,
                internal::RealBox<typename array_traits<Var_t>::Scalar_t,
                                  array_traits<Var_t>::sizeCT, array_traits<Var_t>::containerType,
                                  BS, isFixedRange, MinCT, MaxCT, LRCT>,
                _iFun_, _fFun_, internal::DefaultElectron<Var_t, double>, fOpt, rOpt> {
  using Base_t = typename internal::AOSBase<
      Var_t, double, Arg_t,
      internal::RealBox<typename array_traits<Var_t>::Scalar_t, array_traits<Var_t>::sizeCT,
                        array_traits<Var_t>::containerType, BS, isFixedRange, MinCT, MaxCT, LRCT>,
      _iFun_, _fFun_, internal::DefaultElectron<Var_t, double>, fOpt, rOpt>;

 public:
  HEU_MAKE_AOSBOXED_TYPES(Base_t);

  friend class internal::AOSBase<
      Var_t, double, Arg_t,
      internal::RealBox<typename array_traits<Var_t>::Scalar_t, array_traits<Var_t>::sizeCT,
                        array_traits<Var_t>::containerType, BS, isFixedRange, MinCT, MaxCT, LRCT>,
      _iFun_, _fFun_, internal::DefaultElectron<Var_t, double>, fOpt,
      RecordOption::DONT_RECORD_FITNESS>;

  HEU_RELOAD_MEMBERFUCTION_RUN

 protected:
  static inline bool electronIteratorCompareFun(const ElectronIt_t& a,
                                                const ElectronIt_t& b) noexcept {
    if constexpr (fOpt == FitnessOption::FITNESS_GREATER_BETTER) {
      return a->energy > b->energy;
    } else {
      return a->energy < b->energy;
    }
  }

  void __impl_selectAndMakeLayers() {
    std::vector<ElectronIt_t> elecSortSpace;
    elecSortSpace.reserve(this->_electrons.size());

    for (auto it = this->_electrons.begin(); it != this->_electrons.end(); ++it) {
      elecSortSpace.emplace_back(it);
    }

    std::sort(elecSortSpace.begin(), elecSortSpace.end(), electronIteratorCompareFun);

    // remove worse electrons
    while (elecSortSpace.size() > this->_option.electronNum) {
      this->_electrons.erase(elecSortSpace.back());
      elecSortSpace.pop_back();
    }

    this->_layers.clear();

    // a estimation value of layer number
    const int tempLayerNum = randIdx(1, int(this->_option.maxLayerNum + 1));

    // the accumulate distribution of electrons
    Eigen::ArrayXi numOfEachLayer(tempLayerNum);
    {
      Eigen::ArrayXd floatNumOfEachLayer(tempLayerNum);
      for (int idx = 0; idx < tempLayerNum; idx++) {
        floatNumOfEachLayer[idx] = gaussianCurve<double>(idx + 1, 0, tempLayerNum / 6.0);
      }
      floatNumOfEachLayer /= floatNumOfEachLayer.sum();
      floatNumOfEachLayer *= this->_electrons.size();
      numOfEachLayer = floatNumOfEachLayer.round().cast<int>();
    }

    for (int idx = 1; idx < numOfEachLayer.size(); idx++) {
      numOfEachLayer[idx] += numOfEachLayer[idx - 1];
    }

    for (int& val : numOfEachLayer) {
      if (val > this->_electrons.size()) {
        val = this->_electrons.size();
      }
    }

    this->_layers.resize(tempLayerNum);

    for (Layer_t& layer : this->_layers) {
      layer.clear();
    }

    int curLayerIdx = 0;
    this->_layers.front().reserve(numOfEachLayer[0]);
    for (int countedElectrons = 0; countedElectrons < this->_electrons.size(); countedElectrons++) {
      if (countedElectrons >= numOfEachLayer[curLayerIdx]) {
        curLayerIdx++;
        this->_layers[curLayerIdx].reserve(numOfEachLayer[curLayerIdx] -
                                           numOfEachLayer[curLayerIdx - 1]);
      }

      this->_layers[curLayerIdx].emplace_back(&*elecSortSpace[countedElectrons]);
    }

    // remove empty layers
    while (this->_layers.back().size() <= 0) {
      this->_layers.pop_back();
    }
  }
};

}  // namespace heu

#endif  //  HEU_AOS_HPP