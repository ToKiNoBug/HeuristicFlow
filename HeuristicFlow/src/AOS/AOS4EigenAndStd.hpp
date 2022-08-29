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

#ifndef HEU_AOS4EIGENANDSTD_HPP
#define HEU_AOS4EIGENANDSTD_HPP

#include "InternalHeaderCheck.h"
#include "AOSBoxed.hpp"

namespace heu {
namespace internal {

template <typename Var_t, class Fitness_t, class Arg_t, class Box_t,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::iFun_t _iFun_,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::fFun_t _fFun_, class Electron,
          FitnessOption fOpt>
class AOS4Eigen : public AOSBoxed<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron> {
  using Base_t = AOSBoxed<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron>;

 public:
  HEU_MAKE_AOSBOXED_TYPES(Base_t);

 protected:
  void __impl_computeAtomBSBELE() noexcept {
    auto prevBestPtr = this->_atomBestPtr;
    this->_bindingState.setZero(this->_electrons.front().state.rows(),
                                this->_electrons.front().state.cols());
    this->_bindingEnergy = 0;
    this->_atomBestPtr = &this->_electrons.front();

    for (const Electron_t& elec : this->_electrons) {
      this->_bindingState += elec.state;
      this->_bindingEnergy += elec.energy;
      if constexpr (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
        if (elec.energy < this->_atomBestPtr->energy) {
          this->_atomBestPtr = &elec;
        }
      } else {
        if (elec.energy > this->_atomBestPtr->energy) {
          this->_atomBestPtr = &elec;
        }
      }
    }

    this->_bindingState /= double(this->_electrons.size());
    this->_bindingEnergy /= this->_electrons.size();

    if (prevBestPtr->energy == this->_atomBestPtr->energy) {
      this->_earlyStopCounter++;
    } else {
      this->_earlyStopCounter = 0;
    }
  }

  void __impl_computeLayerBSBELE() noexcept {
    for (int layerIdx = 0; layerIdx < this->_layers.size(); layerIdx++) {
      Layer_t& layer = this->_layers[layerIdx];
      layer.bindingState = layer.front()->state;
      layer.bindingEnergy = layer.front()->energy;
      layer.layerBestIdx = 0;
      for (int idx = 1; idx < layer.size(); idx++) {
        layer.bindingState += layer.at(idx)->state;
        layer.bindingEnergy += layer.at(idx)->energy;
        if constexpr (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
          if (layer.at(idx) < layer.at(layer.layerBestIdx)) {
            layer.layerBestIdx = idx;
          }
        } else {
          if (layer.at(idx) > layer.at(layer.layerBestIdx)) {
            layer.layerBestIdx = idx;
          }
        }
      }

      layer.bindingState /= double(layer.size());
      layer.bindingEnergy /= layer.size();
    }
  }

  void __impl2_applyPhotonEffect(const Electron_t& parent, const Layer_t& layer, const int layerIdx,
                                 Electron_t* child) const noexcept {
    const Var_t& parentState = parent.state;
    Var_t alpha(parentState.rows(), parentState.cols()),
        beta(parentState.rows(), parentState.cols()), gamma(parentState.rows(), parentState.cols());

    randD(alpha.data(), int(alpha.size()));
    randD(beta.data(), int(beta.size()));
    randD(gamma.data(), int(gamma.size()));

    if (Base_t::template isBetter<fOpt>(parent.energy, layer.bindingEnergy)) {
      // E_i^k<BE^k
      child->state =
          parent.state + alpha * (beta * layer.layerBestState() - gamma * layer.bindingState);
    } else {
      // E_i^k>=BE^k
      child->state = parent.state +
                     alpha * (beta * this->bestElectron().state - gamma * this->bindingState()) /
                         (layerIdx + 1);
    }

    for (int varIdx = 0; varIdx < alpha.size(); varIdx++) {
      this->applyConstraint(&child->state, varIdx);
    }
  }

  inline void __impl2_applyNonPhotonEffect(const Electron& parent,
                                           Electron_t* child) const noexcept {
    child->state =
        parent.state + this->delta() * Var_t::Random(parent.state.rows(), parent.state.cols());

    for (int idx = 0; idx < parent.state.size(); idx++) {
      this->applyConstraint(&child->state, idx);
    }
  }
};

template <typename Var_t, class Fitness_t, class Arg_t, class Box_t,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::iFun_t _iFun_,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::fFun_t _fFun_, class Electron,
          FitnessOption fOpt>
class AOS4Std : public AOSBoxed<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron> {
  using Base_t = AOSBoxed<Var_t, Fitness_t, Arg_t, Box_t, _iFun_, _fFun_, Electron>;

 public:
  HEU_MAKE_AOSBOXED_TYPES(Base_t);

 protected:
  void __impl_computeAtomBSBELE() noexcept {
    auto prevBestPtr = this->_atomBestPtr;
    this->_bindingState = this->_electrons.front().state;
    for (auto& val : this->_bindingState) {
      val = 0;
    }

    this->_bindingEnergy = 0;
    this->_atomBestPtr = &this->_electrons.front();

    for (const Electron_t& elec : this->_electrons) {
      // this->_bindingState += elec.state;
      for (int varIdx = 0; varIdx < this->_bindingState.size(); varIdx++) {
        this->_bindingState.operator[](varIdx) += elec.state[varIdx];
      }

      this->_bindingEnergy += elec.energy;
      if constexpr (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
        if (elec.energy < this->_atomBestPtr->energy) {
          this->_atomBestPtr = &elec;
        }
      } else {
        if (elec.energy > this->_atomBestPtr->energy) {
          this->_atomBestPtr = &elec;
        }
      }
    }

    // this->_bindingState /= electrons.size();
    for (auto& val : this->_bindingState) {
      val /= this->_electrons.size();
    }
    this->_bindingEnergy /= this->_electrons.size();

    if (prevBestPtr->energy == this->_atomBestPtr->energy) {
      this->_earlyStopCounter++;
    } else {
      this->_earlyStopCounter = 0;
    }
  }

  void __impl_computeLayerBSBELE() noexcept {
    for (int layerIdx = 0; layerIdx < this->_layers.size(); layerIdx++) {
      Layer_t& layer = this->_layers[layerIdx];
      layer.bindingState = layer.front()->state;
      layer.bindingEnergy = layer.front()->energy;
      layer.layerBestIdx = 0;

      const int varDim = layer.front()->state.size();

      for (int idx = 1; idx < layer.size(); idx++) {
        // layer.bindingState += layer.at(idx)->state;
        for (int varIdx = 0; varIdx < varDim; varIdx++) {
          layer.bindingState[varIdx] += layer.at(idx)->state[varIdx];
        }

        if constexpr (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
          if (layer.at(idx)->energy < layer.layerBestEnergy()) {
            layer.layerBestIdx = idx;
          }
        } else {
          if (layer.at(idx)->energy > layer.layerBestEnergy()) {
            layer.layerBestIdx = idx;
          }
        }
      }

      // layer.bindingState /= layer.size();
      for (auto& var : layer.bindingState) {
        var /= layer.size();
      }
      layer.bindingEnergy /= layer.size();
    }
  }

  void __impl2_applyPhotonEffect(const Electron_t& parent, const Layer_t& layer, const int layerIdx,
                                 Electron_t* child) const noexcept {
    const Var_t& parentState = parent.state;

    Var_t alpha, beta, gamma;
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      alpha.resize(this->dimensions());
      beta.resize(this->dimensions());
      gamma.resize(this->dimensions());
      child->state.resize(this->dimensions());
    }

    randD(alpha.data(), alpha.size());
    randD(beta.data(), beta.size());
    randD(gamma.data(), gamma.size());

    if (Base_t::template isBetter<fOpt>(parent.energy, layer.bindingEnergy)) {
      // E_i^k<BE^k
      for (int varIdx = 0; varIdx < this->dimensions(); varIdx++) {
        child->state[varIdx] =
            parent.state[varIdx] + alpha[varIdx] * (beta[varIdx] * layer.layerBestState()[varIdx] -
                                                    gamma[varIdx] * layer.bindingState[varIdx]);
      }
    } else {
      // E_i^k>=BE^k
      for (int varIdx = 0; varIdx < this->dimensions(); varIdx++) {
        child->state[varIdx] =
            parent.state[varIdx] + alpha[varIdx] *
                                       (beta[varIdx] * this->bestElectron().state[varIdx] -
                                        gamma[varIdx] * this->bindingState()[varIdx]) /
                                       (layerIdx + 1);
      }
    }

    for (int idx = 0; idx < this->dimensions(); idx++) {
      this->applyConstraint(&child->state, idx);
    }
  }

  inline void __impl2_applyNonPhotonEffect(const Electron& parent,
                                           Electron_t* child) const noexcept {
    /*
    child->state =
        parent + this->delta() * Var_t::Random(parent.state.rows(), parent.state.cols());
        */
    child->state = parent.state;
    for (int idx = 0; idx < this->dimensions(); idx++) {
      child->state[idx] += randD(-this->delta(idx), this->delta(idx));

      this->applyConstraint(&child->state, idx);
    }
  }
};
}  // namespace internal
}  // namespace heu

#endif  //  HEU_AOS4EIGENANDSTD_HPP