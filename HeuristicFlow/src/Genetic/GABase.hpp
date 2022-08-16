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

#ifndef HEU_GABASE_H
#define HEU_GABASE_H

#include <vector>
#include <list>
#include <cmath>
#include <random>
#include <algorithm>
#include <type_traits>

#ifdef HEU_DO_OUTPUT
#include <iostream>
#endif

#include <HeuristicFlow/Global>
#include "InternalHeaderCheck.h"
#include "GAOption.hpp"
#include "GAAbstract.hpp"

#include "IsGene.hpp"

namespace heu {

namespace internal {

#if __cplusplus >= 202002L
template <class Gene_t>
concept isGAGene = requires(Gene_t *g, const Gene_t *cg) {
  {is_GA_gene_v<Gene_t>};
  {g->set_fitness_uncomputed()};
};
#endif  //  #if __cplusplus >= 202002L

/**
 * \ingroup HEU_GENETIC
 * \class GABase
 * \brief Genetic algorithm base class.
 *  It's an abstrcat base class for all genetic algorithm solvers.
 *
 * This class maintance its GAOption as a member.It also implements
 *
 * \tparam Var_t Encoded solution type. For instance, when solving TSP problems, what we really want
 * to solve is a permulation to arrange order between nodes, but the permulation is encoded into a
 * double array and decoded via sorting. So in this instance, Var_t is assigned to
 * std::vector<double> instead of an permulation type. You should decode when calculating fitness
 * value. \tparam Fitness_t Type of fitness \tparam Record Whether the solver records the changes of
 * fitness value or not. \tparam Args_t Extra custom parameters. \tparam _iFun_ Function to
 * initialize an individual in population. This function will be called only when initializing. Use
 * nullptr if you hope to determine it at runtime \tparam _fFun_ Funtion to compute fitness for any
 * individual. Use nullptr if you hope to determine it at runtime \tparam _cFun_ Function to apply
 * crossover. Use nullptr if you hope to determine it at runtime \tparam _mFun_ Function to apply
 * mutation. Use nullptr if you hope to determine it at runtime
 */
template <typename Var_t, typename Fitness_t, RecordOption Record, class Gene, class Args_t,
          typename GAAbstract<Var_t, Fitness_t, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Fitness_t, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Fitness_t, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Fitness_t, Args_t>::mutateFun _mFun_>
#if __cplusplus >= 202002L
requires isGAGene<Gene>
#endif  //  #if __cplusplus >= 202002L
class GABase : public GAAbstract<Var_t, Fitness_t, Args_t>,
               public GAAbstract<Var_t, Fitness_t, Args_t>::template iFunBody<_iFun_>,
               public GAAbstract<Var_t, Fitness_t, Args_t>::template fFunBody<_fFun_>,
               public GAAbstract<Var_t, Fitness_t, Args_t>::template cFunBody<_cFun_>,
               public GAAbstract<Var_t, Fitness_t, Args_t>::template mFunBody<_mFun_> {
 private:
  using Base_t = GAAbstract<Var_t, Fitness_t, Args_t>;

 public:
  HEU_MAKE_GAABSTRACT_TYPES(Base_t)

  ~GABase() = default;

 protected:
  using poplist_t = std::list<Gene>;

 public:
  /// Type of gene
  using Gene_t = Gene;
  /// list iterator to Gene
  using GeneIt_t = typename poplist_t::iterator;
  /**
   * \brief Set the option object
   *
   * \param o GAOption
   */
  inline void setOption(const GAOption &o) noexcept { _option = o; }

  /**
   * \brief Get the option object
   *
   * \return const GAOption& Const reference to _option
   */
  inline const GAOption &option() const noexcept { return _option; }

  /**
   * \brief Initialize the whole population
   *
   */
  void initializePop() noexcept {
    _population.resize(_option.populationSize);
    for (auto &i : _population) {
      GAExecutor<Base_t::HasParameters>::doInitialization(this, &i.decision_variable);

      i.set_fitness_uncomputed();
    }
  }

  /**
   * \brief Get the whole population
   *
   * \return const std::list<Gene>& Const reference the the population
   */
  inline const poplist_t &population() const noexcept { return _population; }

  /**
   * \brief Get the population that has been used.
   *
   * \return size_t generation
   */
  inline size_t generation() const noexcept { return _generation; }

  /**
   * \brief Get the current fail times
   *
   * \return size_t fail times
   */
  inline size_t failTimes() const noexcept { return _failTimes; }

 protected:
  poplist_t _population;  ///< Population stored in list
  GAOption _option;       ///< Option of GA solver

  size_t _generation;  ///< Current generation
  size_t _failTimes;   ///< Current failtimes

  inline void __impl_clearRecord() noexcept {
  }  ///< Nothing is need to do if the solver doesn't record fitnesses.

  template <class this_t>
  inline void __impl_recordFitness() noexcept {
  }  ///< Nothing is need to do if the solver doesn't record fitnesses.

  /**
   * \brief Run the solver
   *
   * \tparam this_t Type of a solver. This type is added to record fitness through some template
   * tricks like CRTP.
   */
  template <class this_t>
  void __impl_run() noexcept {
    _generation = 0;
    _failTimes = 0;
    static_cast<this_t *>(this)->__impl_clearRecord();
    while (true) {
      _generation++;
      static_cast<this_t *>(this)->__impl_computeAllFitness();

      static_cast<this_t *>(this)->__impl_select();

      static_cast<this_t *>(this)->template __impl_recordFitness<this_t>();

      if (_generation > _option.maxGenerations) {
#ifdef HEU_DO_OUTPUT
        std::cout << "Terminated by max generation limitation" << std::endl;
#endif
        break;
      }
      if (_option.maxFailTimes > 0 && _failTimes > _option.maxFailTimes) {
#ifdef HEU_DO_OUTPUT
        std::cout << "Terminated by max failTime limitation" << std::endl;
#endif
        break;
      }
#ifdef HEU_DO_OUTPUT
      std::cout << "Generation "
                << _generation
                //<<" , elite fitness="<<_eliteIt->fitness()
                << std::endl;
#endif
      static_cast<this_t *>(this)->__impl_crossover();
      static_cast<this_t *>(this)->__impl_mutate();
    }
    _generation--;
  }

  /**
   * \brief Compute fitness for the whole population
   *
   * If OpenMP is used, this process will be parallelized.
   *
   */
  void __impl_computeAllFitness() noexcept {
#ifdef HEU_HAS_OPENMP
    std::vector<Gene *> tasks;
    tasks.resize(0);
    tasks.reserve(_population.size());
    for (Gene &i : _population) {
      if (i.is_fitness_computed) {
        continue;
      }
      tasks.emplace_back(&i);
    }
    static const int32_t thN = threadNum();
#pragma omp parallel for schedule(dynamic, tasks.size() / thN)
    for (int i = 0; i < int(tasks.size()); i++) {
      Gene *ptr = tasks[i];

      GAExecutor<Base_t::HasParameters>::doFitness(this, &ptr->decision_variable, &ptr->fitness);

      ptr->is_fitness_computed = true;
    }
#else
    for (Gene &i : _population) {
      if (i.is_fitness_computed) {
        continue;
      }

      GAExecutor<Base_t::HasParameters>::doFitness(this, &i.decision_variable, &i.fitness);

      i.is_fitness_computed = true;
    }
#endif
  }

  /**
   * \brief Virtual function to apply selection.
   *
   * Since selection differes a lot, this function can only be implemented by derived classes.
   *
   * Usually the population will exceeds the population size, so this procedure erases relatively
   * worse genes to ensure the size of population can be restored.
   */
  // void select() = 0;

  /**
   * \brief Apply crossover.
   *
   * Apply crossover operation which let randomly 2 gene born 2 more new one. They will be added to
   * population.
   */
  void __impl_crossover() noexcept {
    std::vector<GeneIt_t> crossoverQueue;
    crossoverQueue.clear();
    crossoverQueue.reserve(_population.size());

    for (GeneIt_t it = _population.begin(); it != _population.end(); it++) {
      if (randD() <= _option.crossoverProb) {
        crossoverQueue.emplace_back(it);
      }
    }

    std::shuffle(crossoverQueue.begin(), crossoverQueue.end(), global_mt19937());

    if (crossoverQueue.size() % 2 == 1) {
      crossoverQueue.pop_back();
    }

    while (!crossoverQueue.empty()) {
      GeneIt_t a, b;
      a = crossoverQueue.back();
      crossoverQueue.pop_back();
      b = crossoverQueue.back();
      crossoverQueue.pop_back();
      _population.emplace_back();
      Gene *childA = &_population.back();
      _population.emplace_back();
      Gene *childB = &_population.back();

      GAExecutor<Base_t::HasParameters>::doCrossover(
          this, &a->decision_variable, &b->decision_variable, &childA->decision_variable,
          &childB->decision_variable);

      childA->set_fitness_uncomputed();
      childB->set_fitness_uncomputed();
    }
  }

  /**
   * \brief Apply mutation operation which slightly modify gene. The modified gene will add to the
   * population.
   */
  void __impl_mutate() noexcept {
    std::vector<GeneIt_t> mutateList;
    mutateList.reserve(size_t(this->_population.size() * this->_option.mutateProb * 2));
    for (auto it = this->_population.begin(); it != this->_population.end(); ++it) {
      if (randD() <= this->_option.mutateProb) {
        mutateList.emplace_back(it);
      }
    }
    for (auto src : mutateList) {
      this->_population.emplace_back();
      GAExecutor<Base_t::HasParameters>::doMutation(this, &src->decision_variable,
                                                    &this->_population.back().decision_variable);
      this->_population.back().set_fitness_uncomputed();
    }
  }

 protected:
  // Internal class to apply the four operation
  template <bool HasParameters, class unused = void>
  struct GAExecutor {
    inline static void doInitialization(GABase *s, Var_t *v) noexcept { s->runiFun(v, &s->_args); }

    inline static void doFitness(GABase *s, const Var_t *v, Fitness_t *f) noexcept {
      s->runfFun(v, &s->_args, f);
    }

    inline static void doCrossover(GABase *s, const Var_t *p1, const Var_t *p2, Var_t *c1,
                                   Var_t *c2) noexcept {
      s->runcFun(p1, p2, c1, c2, &s->_args);
    }

    inline static void doMutation(GABase *s, const Var_t *src, Var_t *dst) noexcept {
      s->runmFun(src, dst, &s->_args);
    }
    static_assert(HasParameters == GABase::HasParameters,
                  "struct GAExecutor actived with wrong template parameter.");
  };

  // Internal class to apply the four operation
  template <class unused>
  struct GAExecutor<false, unused> {
    inline static void doInitialization(GABase *s, Var_t *v) noexcept { s->runiFun(v); }

    inline static void doFitness(GABase *s, const Var_t *v, Fitness_t *f) noexcept {
      s->runfFun(v, f);
    }

    inline static void doCrossover(GABase *s, const Var_t *p1, const Var_t *p2, Var_t *c1,
                                   Var_t *c2) noexcept {
      s->runcFun(p1, p2, c1, c2);
    }

    inline static void doMutation(GABase *s, const Var_t *src, Var_t *dst) noexcept {
      s->runmFun(src, dst);
    }

    static_assert(GABase::HasParameters == false,
                  "struct GAExecutor actived with wrong template parameter.");
  };
};

#define HEU_MAKE_GABASE_TYPES(Base_t)         \
  using Gene_t = typename Base_t::Gene_t;     \
  using GeneIt_t = typename Base_t::GeneIt_t; \
  HEU_MAKE_GAABSTRACT_TYPES(Base_t)

/**
 * \ingroup HEU_GENETIC
 * \class GABase
 * \brief partial specialization for GABase with record.
 *
 * GABase with fitness maintains a record of fitness.
 * It's an abstrcat base class for all genetic algorithm solvers that records fitness.
 *
 * \tparam Var_t  Type of decisition variable
 * \tparam Fitness_t  Type of fitness value(objective value)
 * \tparam RecordOption  Whether the solver records fitness changelog
 * \tparam Args_t  Type of other parameters.
 *
 * \note GABase with record is derived from GABase without record to avoid duplicated
 * implementation.
 */
template <typename Var_t, typename Fitness_t, class Gene, class Args_t,
          typename GAAbstract<Var_t, Fitness_t, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Fitness_t, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Fitness_t, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Fitness_t, Args_t>::mutateFun _mFun_>
class GABase<Var_t, Fitness_t, RECORD_FITNESS, Gene, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>
    : public GABase<Var_t, Fitness_t, DONT_RECORD_FITNESS, Gene, Args_t, _iFun_, _fFun_, _cFun_,
                    _mFun_> {
 private:
  using Base_t =
      GABase<Var_t, Fitness_t, DONT_RECORD_FITNESS, Gene, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;
  friend class GABase<Var_t, Fitness_t, DONT_RECORD_FITNESS, Gene, Args_t, _iFun_, _fFun_, _cFun_,
                      _mFun_>;

 public:
  HEU_MAKE_GABASE_TYPES(Base_t)

  ~GABase() noexcept {}

  /// best fitness

  /// result

  /**
   * \brief Get the record.
   *
   * \return const std::vector<Fitness_t> & Const reference to the record.
   */
  const std::vector<Fitness_t> &record() const noexcept { return _record; }

 protected:
  /**
   * \brief Record the best fitness of each generation.
   *
   * \note This member only exists when RecordOption is RECORD_FITNESS
   */
  std::vector<Fitness_t> _record;

  inline void __impl_clearRecord() noexcept {
    _record.clear();
    _record.reserve(this->_option.maxGenerations + 1);
  }

  template <class this_t>
  inline void __impl_recordFitness() noexcept {
    _record.emplace_back(static_cast<this_t *>(this)->bestFitness());
  }
};
}  //  namespace internal
}  //  namespace heu

#endif  // HEU_GABASE_H
