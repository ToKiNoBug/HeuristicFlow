#ifndef HEU_AOSBOXED_HPP
#define HEU_AOSBOXED_HPP

#include <list>
#include <HeuristicFlow/SimpleMatrix>

#include "AOSParameterPack.hpp"
#include "AOSOption.hpp"

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

  using FastFitness_t = typename std::conditional<sizeof(Fitness_t) <= sizeof(void*), Fitness_t,
                                                  const Fitness_t&>::type;

  class Layer : public std::vector<Electron*> {
   public:
    Var_t bindingState;
    Fitness_t bindingEnergy;
    int layerBestIdx;

    inline Var_t& layerBestState() noexcept { return at(layerBestIdx)->state; }

    inline const Var_t& layerBestState() const noexcept { return at(layerBestIdx)->state; }

    inline FastFitness_t layerBestEnergy() const noexcept { return at(layerBestIdx)->energy; }
  };

  inline AOSOption& option() noexcept { return _option; }

  inline const AOSOption& option() const noexcept { return _option; }

  inline void setOption(const AOSOption& _opt) noexcept { _option = _opt; }

  inline const Electron_t& bestElectron() const { return *_atomBestPtr; }

  inline const Var_t& bindingState() const { return _bindingState; }

  inline FastFitness_t bindingEnergy() const { return _bindingEnergy; }

  inline const std::list<Electron_t>& electrons() const {return _electrons};

  inline const std::vector<Layer>& layers() const {return _layers};

 protected:
  std::list<Electron_t> _electrons;
  std::vector<Layer> _layers;
  Var_t _bindingState;
  Fitness_t _bindingEnergy;
  const Electron_t* _atomBestPtr;
  AOSOption _option;

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

  template <Fitness fOpt>
  struct LayerExecutor<false, fOpt> {
    static inline void updateLayerBSBELE(layer* layer) {}
  };

  template <FitnessOpt fOpt>
  inline static bool isBetter(FastFitness_t a, FastFitness_t b) noexcept {
    if (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
      return a < b;
    }
    return a > b;
  }
};

#define HEU_MAKE_AOSBOXED_TYPES(Base_t)         \
  HEU_MAKE_AOSPARAMETERPACK_TYPES(Base_t)       \
  using Electron_t = typename Base_t::Electron; \
  using Layer_t = typename Base_t::Layer;       \
  using FastFitness_t = typename Base_t ::FastFitness_t;

}  // namespace internal
}  // namespace heu

#endif  //  HEU_AOSBOXED_HPP