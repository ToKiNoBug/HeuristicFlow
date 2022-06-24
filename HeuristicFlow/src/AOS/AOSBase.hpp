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

#ifndef HEU_AOSBASE_HPP
#define HEU_AOSBASE_HPP

#include "InternalHeaderCheck.hpp"
#include "AOS4EigenAndStd.hpp"

namespace heu {
namespace internal {

template <typename Var_t, class Fitness_t, class Arg_t, class Box_t,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::iFun_t _iFun_,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::fFun_t _fFun_, class Electron,
          FitnessOption fOpt, RecordOption rOpt>
class AOSBase
    : public std::conditional<
          array_traits<Var_t>::isEigenClass,
          class AOS4Eigen<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron, fOpt>,
          class AOS4Std<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron, fOpt>>::type {
  using Base_t = typename std::conditional<
      array_traits<Var_t>::isEigenClass,
      class AOS4Eigen<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron, fOpt>,
      class AOS4Std<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron, fOpt>>::type;

 public:
  HEU_MAKE_AOSBOXED_TYPES(Base_t);

  void initializePop() {
    this->_layers.clear();

    {
      const bool assignedMaxLayerNumMustNotBeGreaterThanElectronNum =
          (this->_option.maxLayerNum <= this->_option.electronNum);
      assert(assignedMaxLayerNumMustNotBeGreaterThanElectronNum);
    }

    this->_electrons.clear();

    this->_electrons.resize(this->_option.electronNum);

    for (Electron_t& elec : this->_electrons) {
      Base_t::template AOSExecutor<>::doInitiailization(this, &elec.state);
      elec.setUncomputed();
    }

    this->_atomBestPtr = &this->_electrons.front();

    this->_generation = 0;
    this->_earlyStopCounter = 0;

    this->_layers.clear();
    this->_layers.resize(this->_option.electronNum);
    for (Layer_t& layer : this->_layers) {
      layer.reserve(this->_option.electronNum);
    }
    this->_layers.clear();
  }

 protected:
  inline void __impl_recordFitness() const noexcept {}

  inline bool __impl_shouldTerminate() const noexcept {
    return this->_generation >= this->_option.maxGeneration ||
           this->_earlyStopCounter >= this->_option.maxEarlyStop;
  }

  template <class this_t>
  inline void __impl_updateElectrons() {
    for (int layerIdx = 0; layerIdx < this->_layers.size(); layerIdx++) {
      Layer_t& curLayer = this->_layers[layerIdx];
      for (Electron_t* elecPtr : this->_layers[layerIdx]) {
        this->_electrons.emplace_back();
        Electron_t* newPtr = &this->_electrons.back();
        if (randD() >= this->_option.photonRate) {
          static_cast<this_t*>(this)->__impl2_applyPhotonEffect(*elecPtr, curLayer, layerIdx,
                                                                newPtr);
        } else {
          static_cast<this_t*>(this)->__impl2_applyNonPhotonEffect(*elecPtr, newPtr);
        }

        newPtr->setUncomputed();
      }
    }
  }

  template <class this_t>
  void __impl_run() {
    while (true) {
      static_cast<this_t*>(this)->__impl_computeFitness();

      static_cast<this_t*>(this)->__impl_computeAtomBSBELE();

      static_cast<this_t*>(this)->__impl_recordFitness();

      static_cast<this_t*>(this)->__impl_selectAndMakeLayers();

      if (static_cast<this_t*>(this)->__impl_shouldTerminate()) {
        break;
      }

      static_cast<this_t*>(this)->__impl_computeLayerBSBELE();

      static_cast<this_t*>(this)->template __impl_updateElectrons<this_t>();

      this->_generation++;
    }
  }

  static_assert(rOpt == RecordOption::DONT_RECORD_FITNESS, "Wrong specilization!");
};

template <typename Var_t, class Fitness_t, class Arg_t, class Box_t,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::iFun_t _iFun_,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::fFun_t _fFun_, class Electron,
          FitnessOption fOpt>
class AOSBase<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron, fOpt,
              RecordOption::RECORD_FITNESS>
    : public AOSBase<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron, fOpt,
                     RecordOption::DONT_RECORD_FITNESS> {
  using Base_t = AOSBase<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron, fOpt,
                         RecordOption::DONT_RECORD_FITNESS>;

 public:
  HEU_MAKE_AOSBOXED_TYPES(Base_t);

  friend Base_t;

  const std::vector<Fitness_t>& record() const noexcept { return _record; }

  inline void initializePop() {
    Base_t::initializePop();
    _record.clear();
    _record.reserve(this->_option.maxGeneration + 1);
  }

 protected:
  std::vector<Fitness_t> _record;

  inline void __impl_recordFitness() { _record.emplace_back(this->_atomBestPtr->energy); }
};

}  // namespace internal
}  // namespace heu

#endif  //  HEU_AOSBASE_HPP