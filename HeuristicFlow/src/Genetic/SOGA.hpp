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

#ifndef HEU_SOGA_H
#define HEU_SOGA_H

#include "InternalHeaderCheck.h"
#include "GABase.hpp"
#include "SOGASelecters.hpp"
#include "DefaultGeneType.hpp"

namespace heu {

/**
 * \ingroup HEU_GENETIC
 * \class SOGA
 * \brief Single-object genetic solver.
 *
 * \tparam Var_t  Type of decisition variable.
 * \tparam fOpt Whether greater fitness value means better. SOGA will always try to find the best,
 * so it's vital to tell SOGA which direction is right. \tparam Record  Whether the solver records
 * fitness changelog. \tparam Args_t  Type of other parameters. \tparam _iFun_ Function to
 * initialize an individual in population. This function will be called only when initializing. Use
 * nullptr if you hope to determine it at runtime. nullptr as default value. \tparam _fFun_ Funtion
 * to compute fitness for any individual. Use nullptr if you hope to determine it at runtime.
 * nullptr as default value. \tparam _cFun_ Function to apply crossover. Use nullptr if you hope to
 * determine it at runtime. nullptr as default value. \tparam _mFun_ Function to apply mutation. Use
 * nullptr if you hope to determine it at runtime. nullptr as default value.
 *
 *
 * \note `Arg_t` is a pseudo-global variable for solvers. It's not global variables because it's
 * stored as a member of solver, thus if in some special conditions you have to run multiple
 * solvers, each solver has its own gloabl variable. It can be any type about the problem you are
 * solving. For constrainted problems, it can be a box-constraint type; for TSP problems, it can be
 * the distance matrix; and even for GA-BP or PSO-BP, it can be your data set. If you don't need it,
 * use `void` and it will disapper.
 *
 * Usually there are five main procedures in GA :
 * 1. initialization : Filling all decision variables in the population with random initial values.
 * 2. computing fitness : Go through the whole population and compute each's fitness value.
 * 3. selection : Eliminate some relatively bad inidividual sothat the size of population can
 * decline to `populationSize` assigned in the `GAOption`.
 * 4. crossover Randomly pick individuals in pairs as parents and create pairs of child
 * inidividuals. The childern will be inserted into the original population.
 * 5. mutation Randomly pick individuals and change their values slightly. This procedure doesn't
 * change the value inplace, but insert the modified value into the population.
 *
 * The 1,2,4,5 each corresponds to a function representively, called initialization function(iFun),
 * fitness function(fFun), crossover function(cFun), mutation function(mFun). These functions can be
 * assigned either at compile time or runtime.
 *
 * APIs of SOGA are implemented sperately in many internal base classes, and they are organized
 * through public inheriting.
 *
 * ## APIs that **all** genetic solvers have:
 * - `void setOption(const GAOption&)` sets the option of solver.
 * - `void initializePop()` that initialized the population.
 * - `void run()` runs the genetic algorithm.
 * - `Fitness_t bestFitness() const` returns the fitness of best solution in current population.
 * - `size_t generation() const` returns the generation that solvers has passed.
 * - `size_t failTimes() const` returns the fail times of current population.
 * - `const std::list<Gene_t> & population() const` returns a const-reference to the population.
 * - `const GAOption & option() const` returns a const-reference to the GAOption of solver.
 * - `typename solver_t::initializeFun` is the type of iFuns
 * - `typename solver_t::fitnessFun` is the type of fFuns
 * - `typename solver_t::crossoverFun` is the type of cFuns
 * - `typename solver_t::mutateFun` is the type of mFuns
 * - `initializeFun iFun() const` returns the initialization function.
 * - `fitnessFun fFun() const` returns the fitness function.
 * - `crossoverFun cFun() const` returns the crossover function.
 * - `mutateFun mFun() const` returns the mutation function.
 * - 'SelectMethod selectMethod() const' returns the selection method.
 *
 *
 * ## APIs that all genetic solvers whose `Args_t` is not `void` have:
 * - `const Arg_t & args() const` returns a const-reference of args.
 * - `void setArgs(const Arg_t &)` set the value of args.
 *
 * ## APIs that all genetic solvers with recording have:
 * - `const std::vector<Fitness_t> & record() const` returns a const reference to the recoding.
 *
 *
 * ## APIs that SOGA solvers have:
 * - `const Var_t& result() const` returns a const-reference the the final result decision variable.
 *
 * \sa GAOption NSGA2 NSGA3 SelectMethod
 */
template <typename Var_t, FitnessOption fOpt = FITNESS_LESS_BETTER,
          RecordOption Record = DONT_RECORD_FITNESS,
          SelectMethod smOpt = SelectMethod::RouletteWheel, class Args_t = void,
          typename internal::GAAbstract<Var_t, double, Args_t>::initializeFun _iFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::fitnessFun _fFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::crossoverFun _cFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::mutateFun _mFun_ = nullptr>
class SOGA : public internal::GABase<Var_t, double, Record, internal::DefaultGene_t<Var_t, double>,
                                     Args_t, _iFun_, _fFun_, _cFun_, _mFun_>,
             public internal::SOGASelector<smOpt> {
 private:
  using Base_t = internal::GABase<Var_t, double, Record, internal::DefaultGene_t<Var_t, double>,
                                  Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;
  friend class internal::GABase<Var_t, double, DONT_RECORD_FITNESS,
                                internal::DefaultGene_t<Var_t, double>, Args_t, _iFun_, _fFun_,
                                _cFun_, _mFun_>;

  template <SelectMethod _sm>
  friend class internal::SOGASelector;

 public:
  HEU_MAKE_GABASE_TYPES(Base_t)
  SOGA() {
    if constexpr (smOpt == SelectMethod::Boltzmann) {
      this->_boltzmannSelectStrength = 10 * (fOpt == FITNESS_LESS_BETTER ? -1 : 1);
    }
  }
  ~SOGA() = default;

  HEU_RELOAD_MEMBERFUCTION_RUN

  static constexpr FitnessOption FitnessOpt = fOpt;

  /**
   * \brief Get the best gene
   *
   * \return const Gene_t& The best gene
   */
  inline const Gene_t& bestGene() const noexcept { return *_bestGene; }

  /**
   * \brief Get the fitness value of best gene
   *
   * \return Fitness_t Best fitness
   */
  inline double bestFitness() const noexcept { return _bestGene->fitness; }

  /**
   * \brief Get result (Var_t).
   *
   * \return const Var_t& The decision variable of the elite gene.
   */
  inline const Var_t& result() const noexcept { return _bestGene->decision_variable; }

  /**
   * \brief Initialize the population and assign the first gene to be the elite.
   *
   */
  inline void initializePop() noexcept {
    Base_t::initializePop();
    _bestGene = this->_population.begin();
  }

 protected:
  GeneIt_t _bestGene;  ///< Iterator the the elite
  /**
   * \brief Returns whether A is better than B
   *
   * \param A The first input
   * \param B The second input
   * \return true A is better than B
   * \return false A is worse than or equals to B
   */
  static inline bool isBetter(double A, double B) noexcept {
    if constexpr (fOpt == FitnessOption::FITNESS_GREATER_BETTER)
      return A > B;
    else
      return A < B;
  }

  static inline bool GeneItCompareFun(const GeneIt_t& a, const GeneIt_t& b) noexcept {
    return isBetter(a->fitness, b->fitness);
  }

  GeneIt_t findCurrentBestGene() noexcept {
    GeneIt_t gIt = this->_population.begin();
    for (GeneIt_t it = this->_population.begin(); it != this->_population.end(); ++it) {
      if (GeneItCompareFun(it, gIt)) {
        gIt = it;
      }
    }

    return gIt;
  }

  inline void applyErasement(const std::list<std::pair<GeneIt_t, double>>& eraseList) noexcept {
    // erase all eliminated candidates from the linked list
    for (auto& pair : eraseList) {
      this->_population.erase(pair.first);
    }
  }

  inline void updateFailTimesAndBestGene(const GeneIt_t& newBestGeneIt,
                                         const double prevFitess) noexcept {
    if (!isBetter(newBestGeneIt->fitness, prevFitess)) {
      this->_failTimes++;
      _bestGene = newBestGeneIt;
    } else {
      this->_failTimes = 0;
      _bestGene = newBestGeneIt;
    }
  }

  /**
   * \brief Call this function only when you are sure that `_bestGene` will always be vaild AFTER
   * the selection.
   */
  inline void updateFailTimesAndBestGene(const GeneIt_t& newBestGeneIt) noexcept {
    updateFailTimesAndBestGene(newBestGeneIt, this->_bestGene->fitness);
  }

  /**
   * \brief Simple select implementation.
   *
   * This function sorts the whole population with their fitness, and numbers of original
   * populationSize genes will be reserved. The rest will be eliminated. The best gene in population
   * will be assigned to be elite.
   *
   */
  inline void __impl_select() noexcept { this->template __impl___impl_select<SOGA>(); }
};

}  //    namespace heu

#endif  // HEU_SOGA_H
