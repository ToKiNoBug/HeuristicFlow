#ifndef HEU_AOSBOXED_HPP
#define HEU_AOSBOXED_HPP

#include "AOSParameterPack.hpp"
#include <list>

namespace heu {
namespace internal {

template <typename Var_t, class Fitness_t, class Arg_t, class Box_t,
          AOSParameterPack<Var_t, Fitness_t, Arg_t>::iFun_t _iFun_,
          AOSParameterPack<Var_t, Fitness_t, Arg_t>::fFun_t _fFun_, class Electron>
class AOSBoxed : public AOSParameterPack<Var_t, Fitness_t, Arg_t>,
                 AOSParameterPack<Var_t, Fitness_t, Arg_t>::template iFunBody<_iFun_>,
                 AOSParameterPack<Var_t, Fitness_t, Arg_t>::template fFunBody<_fFun_>,
                 Box_t {
  using Base_t = AOSParameterPack<Var_t, Fitness_t, Arg_t>;

 public:
  HEU_MAKE_AOSPARAMETERPACK_TYPES(Base_t);
  using Electron_t = Electron;

  class Layer : public std::vector<Electron*> {
   public:
    Var_t bindingState;
    Fitness_t bindingEnergy;
    int layerBestIdx;
  };

 protected:
  std::list<Electron_t> _electrons;
  void __impl_computeFitness() {
#ifdef EIGEN_HAS_OPENMP
    std::vector<Electron_t*> tasks;
    tasks.resize(0);
    tasks.reserve(_electrons.size());
    for (Electron_t& i : _electrons) {
      if (i.isComputed) {
        continue;
      }
      tasks.emplace_back(&i);
    }
    static const int32_t thN = Eigen::nbThreads();
#pragma omp parallel for schedule(dynamic, tasks.size() / thN)
    for (int i = 0; i < tasks.size(); i++) {
      Electron_t* ptr = tasks[i];

      AOSExecutor<Base_t::hasParameters>::doFitness(this, &ptr->self, &ptr->_Fitness);

      ptr->isComputed = true;
    }
#else
    for (Electron_t& i : _electrons) {
      if (i.isComputed) {
        continue;
      }

      AOSExecutor<Base_t::hasParameters>::doFitness(this, &i.self, &i._Fitness);

      i.isComputed = true;
    }
#endif
  }

 protected:
  template <bool hasParameter, typename unused = void>
  struct AOSExecutor {
    static inline void doInitiailization(const AOSBoxed* solver, Var_t* v) {
      solver->runiFun(v, &this->arg());
    }

    static inline void doFitness(const AOSBoxed* solver, const Var_t* v, Fitness* f) {
      solver->runfFun(v, &this->arg(), f);
    }
  };

  template <typename unused>
  struct AOSExecutor<false, unused> {
    static inline void doInitiailization(const AOSBoxed* solver, Var_t* v) { solver->runiFun(v); }

    static inline void doFitness(const AOSBoxed* solver, const Var_t* v, Fitness* f) {
      solver->runfFun(v, f);
    }
  };

  template <bool isEigenClass, FitnessOption fOpt>
  struct LayerExecutor {
    static inline void updateLayerBSBELE(Layer* layer) {
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

    static inline void updateAtomBSBELE(const std::list<Electron_t>& electrons, Var_t* BS,
                                        Fitness_t* BE, const Electron_t** LEIt) {
      BS->setZero(electrons.front().state.rows(), electrons.front().state.cols());
      *BE = 0;
      *LEIt = &electrons.front();

      for (const Electron_t& elec : electrons) {
        *BS += elec.state;
        *BE += elec.energy;
        if (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
          if (elec.energy < *LEIt->energy) {
            *LEIt = &elec;
          }
        } else {
          if (elec.energy > *LEIt->energy) {
            *LEIt = &elec;
          }
        }
      }

      *BS /= electrons.size();
      *BE /= electrons.size();
    }
  };

  template <Fitness fOpt>
  struct LayerExecutor<false, fOpt> {
    static inline void updateLayerBSBELE(layer* layer) {
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

    static inline void updateAtomBSBELE(const std::list<Electron_t>& electrons, Var_t* BS,
                                        Fitness_t* BE, const Electron_t** LEIt) {
      *BS = electrons.front().state;
      for (auto& val : *BS) {
        val = 0;
      }

      *BE = 0;
      *LEIt = &electrons.front();

      for (const Electron_t& elec : electrons) {
        //*BS += elec.state;
        for (int varIdx = 0; varIdx < BS->size(); valIdx++) {
          BS->operator[](valIdx) += elec.energy[valIdx];
        }

        *BE += elec.energy;
        if (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
          if (elec.energy < *LEIt->energy) {
            *LEIt = &elec;
          }
        } else {
          if (elec.energy > *LEIt->energy) {
            *LEIt = &elec;
          }
        }
      }

      //*BS /= electrons.size();
      for (auto& val : *BS) {
        val /= electrons.size();
      }
      *BE /= electrons.size();
    }
  };
};

}  // namespace internal
}  // namespace heu

#endif  //  HEU_AOSBOXED_HPP