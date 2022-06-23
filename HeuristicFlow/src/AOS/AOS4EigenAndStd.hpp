#ifndef HEU_AOS4EIGENANDSTD_HPP
#define HEU_AOS4EIGENANDSTD_HPP

#include "InternalHeaderCheck.hpp"
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
  void __impl_computeAtomBSBELE() {
    this->_bindingState.setZero(this->_electrons.front().state.rows(),
                                this->_electrons.front().state.cols());
    this->_bindingEnergy = 0;
    this->_atomBestPtr = &this->_electrons.front();

    for (const Electron_t& elec : this->_electrons) {
      this->_bindingState += elec.state;
      this->_bindingEnergy += elec.energy;
      if (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
        if (elec.energy < this->_atomBestPtr->energy) {
          this->_atomBestPtr = &elec;
        }
      } else {
        if (elec.energy > this->_atomBestPtr->energy) {
          this->_atomBestPtr = &elec;
        }
      }
    }

    this->_bindingState /= this->_electrons.size();
    this->_bindingEnergy /= this->_electrons.size();
  }

  void __impl_computeLayerBSBELE() {
#ifdef EIGEN_HAS_OPENMP
    static const int32_t thN = Eigen::nbThreads();
#pragma omp parallel for schedule(dynamic, this->_layers.size() / thN)
#endif  //  EIGEN_HAS_OPENMP
    for (int idx = 0; idx < this->_layer.size(); idx++) {
      Layer_t& layer = this->_layers[idx];
      layer->bindingState = layer->front()->state;
      layer->bindingEnergy = layer->front()->energy;
      layer->layerBestIdx = 0;
      for (int idx = 1; idx < layer->size(); idx++) {
        layer->bindingState += layer->at(idx)->state;
        layer->bindingEnergy += layer->at(idx)->energy;
        if (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
          if (layer->at(idx) < layer->at(layer->layerBestIdx)) {
            layer->layerBestIdx = idx;
          }
        } else {
          if (layer->at(idx) > layer->at(layer->layerBestIdx)) {
            layer->layerBestIdx = idx;
          }
        }
      }

      layer->bindingState /= layer->size();
      layer->bindingEnergy /= layer->size();
    }
  }

  void __impl2_applyPhotonEffect(const Electron_t& parent, const Layer_t& layer, const int layerIdx,
                                 Electron_t* child) const {
    const Var_t& parentState = parent.state;
    Var_t alpha(parentState.rows(), parentState.cols()),
        beta(parentState.rows(), parentState.cols()), gamma(parentState.rows(), parentState.cols());

    randD(alpha.data(), alpha.size());
    randD(beta.data(), beta.size());
    randD(gamma.data(), gamma.size());

    if (Base_t::isBetter<fOpt>(parent.energy, layer.bindingEnergy)) {
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
      this->applyConstraint(&child->state(varIdx), varIdx);
    }
  }

  inline void __impl2_applyNonPhotonEffect(const Electron& parent, Electron_t* child) const {
    child->state =
        parent + this->learnRate() * Var_t::Random(parent.state.rows(), parent.state.cols());

    for (int idx = 0; idx < parent.state.size(); idx++) {
      this->applyConstraint(&child->state(idx), idx);
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
  void __impl_computeAtomBSBELE() {
    this->_bindingState = this->_electrons.front().state;
    for (auto& val : this->_bindingState) {
      val = 0;
    }

    this->_bindingEnergy = 0;
    this->_atomBestPtr = &this->_electrons.front();

    for (const Electron_t& elec : this->_electrons) {
      // this->_bindingState += elec.state;
      for (int varIdx = 0; varIdx < this->_bindingState.size(); varIdx++) {
        this->_bindingState.operator[](varIdx) += elec.energy[varIdx];
      }

      this->_bindingEnergy += elec.energy;
      if (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
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
  }

  void __impl_computeLayerBSBELE() {
#ifdef EIGEN_HAS_OPENMP
    static const int32_t thN = Eigen::nbThreads();
#pragma omp parallel for schedule(dynamic, this->_layers.size() / thN)
#endif  //  EIGEN_HAS_OPENMP
    for (int idx = 0; idx < this->_layers.size(); idx++) {
      Layer_t& layer = this->_layers[idx];
      layer->bindingState = layer->front()->state;
      layer->bindingEnergy = layer->front()->energy;
      layer->layerBestIdx = 0;

      const int varDim = layer->front()->size();

      for (int idx = 1; idx < layer->size(); idx++) {
        // layer->bindingState += layer->at(idx)->state;
        for (int varIdx = 0; varIdx < varDim; varIdx++) {
          layer->bindingEnergy[varIdx] += layer->at(idx)->energy[varIdx];
        }

        if (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
          if (layer->at(idx) < layer->at(layer->layerBestIdx)) {
            layer->layerBestIdx = idx;
          }
        } else {
          if (layer->at(idx) > layer->at(layer->layerBestIdx)) {
            layer->layerBestIdx = idx;
          }
        }
      }

      // layer->bindingState /= layer->size();
      for (auto& var : layer->bindingState) {
        var /= layer->size();
      }
      layer->bindingEnergy /= layer->size();
    }
  }

  void __impl2_applyPhotonEffect(const Electron_t& parent, const Layer_t& layer, const int layerIdx,
                                 Electron_t* child) const {
    const Var_t& parentState = parent.state;

#warning Initialization here is waiting to be optimized
    Var_t alpha = parentState, beta = parentState, gamma = parentState;
#warning and here
    child->state = parent.state;

    randD(alpha.data(), alpha.size());
    randD(beta.data(), beta.size());
    randD(gamma.data(), gamma.size());

    if (Base_t::isBetter<fOpt>(parent.energy, layer.bindingEnergy)) {
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
      this->applyConstraint(&child->state[idx], idx);
    }
  }

  inline void __impl2_applyNonPhotonEffect(const Electron& parent, Electron_t* child) const {
    /*
    child->state =
        parent + this->learnRate() * Var_t::Random(parent.state.rows(), parent.state.cols());
        */
    child->state = parent.state;
    for (int idx = 0; idx < this->dimensions(); idx++) {
      child->state[idx] += randD(-this->learningRate(), this->learningRate());
      this->applyConstraint(&child->state[idx], idx);
    }
  }
};
}  // namespace internal
}  // namespace heu

#endif  //  HEU_AOS4EIGENANDSTD_HPP