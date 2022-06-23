#ifndef HEU_AOS_HPP
#define HEU_AOS_HPP

#include <algorithm>

#include "InternalHeaderCheck.hpp"
#include "AOSBase.hpp"
#include "AOSDefaultElectron.hpp"

namespace heu {

template <typename Var_t, BoxShape BS, FitnessOption fOpt = FitnessOption::FITNESS_LESS_BETTER,
          RecordOption rOpt = RecordOption::DONT_RECORD_FITNESS, class Arg_t = void,
          typename internal::AOSParameterPack<Var_t, double, Arg_t>::fFun_t _fFun_ = nullptr,
          bool isFixedRange = false, DivCode MinCT = DivEncode<0, 1>::code,
          DivCode MaxCT = DivEncode<1, 1>::code, DivCode LRCT = DivEncode<1, 50>::code,
          typename internal::AOSParameterPack<Var_t, double, Arg_t>::fFun_t _iFun_ =
              internal::AOSParameterPack<Var_t, double,
                                         Arg_t>::defaultInitializeFunctionThatShouldNotBeCalled>
class AOS : public internal::AOSBase<
                Var_t, double, Arg_t,
                internal::RealBox<array_traits<Var_t>::Scalar_t, array_traits<Var_t>::sizeCT,
                                  array_traits<Var_t>::containerType, BS, isFixedRange, MinCT,
                                  MaxCT, LRCT>,
                _iFun_, _fFun_, internal::DefaultElectron<Var_t, double>, fOpt, rOPt> {
  using Base_t = typename internal::AOSBase<
      Var_t, double, Arg_t,
      internal::RealBox<array_traits<Var_t>::Scalar_t, array_traits<Var_t>::sizeCT,
                        array_traits<Var_t>::containerType, BS, isFixedRange, MinCT, MaxCT, LRCT>,
      _iFun_, _fFun_, internal::DefaultElectron<Var_t, double>, fOpt, rOPt>;

 public:
  HEU_MAKE_AOSBOXED_TYPES(Base_t);

  friend class internal::AOSBase<
      Var_t, double, Arg_t,
      internal::RealBox<array_traits<Var_t>::Scalar_t, array_traits<Var_t>::sizeCT,
                        array_traits<Var_t>::containerType, BS, isFixedRange, MinCT, MaxCT, LRCT>,
      _iFun_, _fFun_, internal::DefaultElectron<Var_t, double>, fOpt,
      RecordOption::DONT_RECORD_FITNESS>;

  HEU_RELOAD_MEMBERFUCTION_RUN

 protected:
  static inline bool electronIteratorCompareFun(const ElectronIt_t& a,
                                                const ElectronIt_t& b) noexcept {
    if (fOpt == FitnessOption::FITNESS_GREATER_BETTER) {
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
    const int tempLayerNum = randIdx(3, int(this->_electrons.size() - 1));

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
    for (int countedElectrons = 0; countedElectrons < this->_electrons.size(); countedElectrons++) {
      if (countedElectrons >= numOfElectrons[curLayerIdx]) {
        curLayerIdx++;
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