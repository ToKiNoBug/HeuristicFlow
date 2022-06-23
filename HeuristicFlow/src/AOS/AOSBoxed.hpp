#ifndef HEU_AOSBOXED_HPP
#define HEU_AOSBOXED_HPP

#include <list>
#include <HeuristicFlow/SimpleMatrix>

#include "InternalHeaderCheck.hpp"
#include "AOSParameterPack.hpp"
#include "AOSOption.hpp"

namespace heu {
namespace internal {

#warning Remember to migrate all functions that compares between Fitness_t s, they are not adpative for potential multi-objective solvers

template <typename Var_t, class Fitness_t, class Arg_t, class Box_t,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::iFun_t _iFun_,
          typename AOSParameterPack<Var_t, Fitness_t, Arg_t>::fFun_t _fFun_, class Electron>
class AOSBoxed : public AOSParameterPack<Var_t, Fitness_t, Arg_t>,
                 AOSParameterPack<Var_t, Fitness_t, Arg_t>::template iFunBody<_iFun_>,
                 AOSParameterPack<Var_t, Fitness_t, Arg_t>::template fFunBody<_fFun_>,
                 Box_t {
  using Base_t = AOSParameterPack<Var_t, Fitness_t, Arg_t>;

 public:
  HEU_MAKE_AOSPARAMETERPACK_TYPES(Base_t);
  using Electron_t = Electron;

  using ElectronIt_t = typename std::list<Electron_t>::iterator;

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

  inline size_t generation() const noexcept { return _generation; }

  inline size_t earlyStopCounter() const noexcept { return _earlyStopCounter; }

 protected:
  std::list<Electron_t> _electrons;
  std::vector<Layer> _layers;
  Var_t _bindingState;
  Fitness_t _bindingEnergy;
  const Electron_t* _atomBestPtr;
  AOSOption _option;
  size_t _generation;
  size_t _earlyStopCounter;

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

      AOSExecutor<>::doFitness(this, &ptr->self, &ptr->_Fitness);

      ptr->isComputed = true;
    }
#else
    for (Electron_t& i : _electrons) {
      if (i.isComputed) {
        continue;
      }

      AOSExecutor<>::doFitness(this, &i.self, &i._Fitness);

      i.isComputed = true;
    }
#endif
  }

 protected:
  template <bool hasParameter = Base_t::hasParameters, typename unused = void>
  struct AOSExecutor {
    static inline void doInitiailization(const AOSBoxed* solver, Var_t* v) {
      if constexpr (Base_t::iFunAtCompileTime ==
                    Base_t::defaultInitializeFunctionThatShouldNotBeCalled) {
        if constexpr (array_traits<Var_t>::isEigenClass) {
          v->resize(solver->dimensions(), 1);
        } else if constexpr (array_traits<Var_t>::sizeCT == -1) {
          v->resize(solver->dimensions());
        }

        for (int idx = 0; idx < solver->dimensions(); idx++) {
          if constexpr (Box_t::Shape == BoxShape::SQUARE_BOX) {
            v->operator[](idx) = randD(this->min(), this->max());
          } else {
            v->operator[](idx) = randD(this->min()[idx], this->max()[idx]);
          }
        }

      } else {
        solver->runiFun(v, &this->arg());
      }
    }

    static inline void doFitness(const AOSBoxed* solver, const Var_t* v, Fitness* f) {
      solver->runfFun(v, &this->arg(), f);
    }
  };

  template <typename unused>
  struct AOSExecutor<false, unused> {
    static inline void doInitiailization(const AOSBoxed* solver, Var_t* v) {
      if constexpr (Base_t::iFunAtCompileTime ==
                    Base_t::defaultInitializeFunctionThatShouldNotBeCalled) {
        if constexpr (array_traits<Var_t>::isEigenClass) {
          v->resize(solver->dimensions(), 1);
        } else if constexpr (array_traits<Var_t>::sizeCT == -1) {
          v->resize(solver->dimensions());
        }

        for (int idx = 0; idx < solver->dimensions(); idx++) {
          if constexpr (Box_t::Shape == BoxShape::SQUARE_BOX) {
            v->operator[](idx) = randD(this->min(), this->max());
          } else {
            v->operator[](idx) = randD(this->min()[idx], this->max()[idx]);
          }
        }

      } else {
        solver->runiFun(v);
      }
    }

    static inline void doFitness(const AOSBoxed* solver, const Var_t* v, Fitness* f) {
      solver->runfFun(v, f);
    }
  };

  template <FitnessOpt fOpt>
  inline static bool isBetter(FastFitness_t a, FastFitness_t b) noexcept {
    if (fOpt == FitnessOption::FITNESS_LESS_BETTER) {
      return a < b;
    }
    return a > b;
  }
};

#define HEU_MAKE_AOSBOXED_TYPES(Base_t)               \
  HEU_MAKE_AOSPARAMETERPACK_TYPES(Base_t)             \
  using Electron_t = typename Base_t::Electron;       \
  using ElectronIt_t = typename Base_t::ElectronIt_t; \
  using Layer_t = typename Base_t::Layer;             \
  using FastFitness_t = typename Base_t ::FastFitness_t;

}  // namespace internal
}  // namespace heu

#endif  //  HEU_AOSBOXED_HPP