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

namespace heu {

/**
 * \ingroup HEU_GENETIC
 * \class SOGA
 * \brief Single-object genetic solver.
 *
 * \tparam Var_t  Type of decisition variable.
 * \tparam fOpt Whether greater fitness value means better. SOGA will always try to find the best, so it's vital to
 * tell SOGA which direction is right.
 * \tparam Record  Whether the solver records fitness changelog.
 * \tparam Args_t  Type of other parameters.
 * \tparam _iFun_ Function to initialize an individual in population. This function will be called only when
 * initializing. Use nullptr if you hope to determine it at runtime. nullptr as default value.
 * \tparam _fFun_ Funtion to compute fitness for any
 * individual. Use nullptr if you hope to determine it at runtime. nullptr as default value.
 * \tparam _cFun_ Function to apply crossover. Use
 * nullptr if you hope to determine it at runtime. nullptr as default value.
 * \tparam _mFun_ Function to apply mutation. Use nullptr if you hope to
 * determine it at runtime. nullptr as default value.
 *
 *
 * \note `Arg_t` is a pseudo-global variable for solvers. It's not global variables because it's stored as a member of
 * solver, thus if in some special conditions you have to run multiple solvers, each solver has its own gloabl
 * variable. It can be any type about the problem you are solving. For constrainted problems, it can be a box-constraint
 * type; for TSP problems, it can be the distance matrix; and even for GA-BP or PSO-BP, it can be your data set. If you
 * don't need it, use `void` and it will disapper.
 *
 * Usually there are five main procedures in GA :
 * 1. initialization : Filling all decision variables in the population with random initial values.
 * 2. computing fitness : Go through the whole population and compute each's fitness value.
 * 3. selection : Eliminate some relatively bad inidividual sothat the size of population can decline to
 * `populationSize` assigned in the `GAOption`.
 * 4. crossover Randomly pick individuals in pairs as parents and create pairs of child inidividuals. The childern will
 * be inserted into the original population.
 * 5. mutation Randomly pick individuals and change their values slightly. This procedure doesn't change the value
 * inplace, but insert the modified value into the population.
 *
 * The 1,2,4,5 each corresponds to a function representively, called initialization function(iFun), fitness
 * function(fFun), crossover function(cFun), mutation function(mFun). These functions can be assigned either at compile
 * time or runtime.
 *
 * APIs of SOGA are implemented sperately in many internal base classes, and they are organized through public
 * inheriting.
 *
 * ## APIs that **all** genetic solvers have:
 * - `void setOption(const GAOption&)` sets the option of solver.
 * - `void initializePop()` that initialized the population.
 * - `void run()` runs the genetic algorithm.
 * - `Fitness_t bestFitness() const` returns the fitness of best solution in current population.
 * - `size_t generation() const` returns the generation that solvers has passed.
 * - `size_t failTimes() const` returns the fail times of current population.
 * - `const std::list<Gene> & population() const` returns a const-reference to the population.
 * - `const GAOption & option() const` returns a const-reference to the GAOption of solver.
 * - `typename solver_t::initializeFun` is the type of iFuns
 * - `typename solver_t::fitnessFun` is the type of fFuns
 * - `typename solver_t::crossoverFun` is the type of cFuns
 * - `typename solver_t::mutateFun` is the type of mFuns
 * - `initializeFun iFun() const` returns the initialization function.
 * - `fitnessFun fFun() const` returns the fitness function.
 * - `crossoverFun cFun() const` returns the crossover function.
 * - `mutateFun mFun() const` returns the mutation function.
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
 * \sa GAOption NSGA2 NSGA3
 */
template <typename Var_t, FitnessOption fOpt = FITNESS_LESS_BETTER, RecordOption Record = DONT_RECORD_FITNESS,
          class Args_t = void, typename internal::GAAbstract<Var_t, double, Args_t>::initializeFun _iFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::fitnessFun _fFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::crossoverFun _cFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::mutateFun _mFun_ = nullptr>
class SOGA : public internal::GABase<Var_t, double, Record, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
 private:
  using Base_t = internal::GABase<Var_t, double, Record, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;
  friend class internal::GABase<Var_t, double, DONT_RECORD_FITNESS, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  HEU_MAKE_GABASE_TYPES(Base_t)
  SOGA() {}
  ~SOGA() {}

  /**
   * \brief Get the fitness value of best gene
   *
   * \return Fitness_t Best fitness
   */
  double bestFitness() const { return _eliteIt->_Fitness; }

  /**
   * \brief Get result (Var_t).
   *
   * \return const Var_t& The decision variable of the elite gene.
   */
  inline const Var_t& result() const { return _eliteIt->self; }

  /**
   * \brief Initialize the population and assign the first gene to be the elite.
   *
   */
  void initializePop() {
    Base_t::initializePop();
    this->_eliteIt = this->_population.begin();
  }

  /**
   * \brief Run the solver.
   *
   */
  inline void run() { this->template __impl_run<SOGA>(); }

 protected:
  GeneIt_t _eliteIt;  ///< Iterator the the elite

  /**
   * \brief Returns whether A is better than B
   *
   * \param A The first input
   * \param B The second input
   * \return true A is better than B
   * \return false A is worse than or equals to B
   */
  static inline bool isBetter(double A, double B) {
    if (fOpt == FitnessOption::FITNESS_GREATER_BETTER) {
      return A > B;
    }
    return A < B;
  }

  /**
   * \brief Simple select implementation.
   *
   * This function sorts the whole population with their fitness, and numbers of original populationSize genes will be
   * reserved. The rest will be eliminated.
   * The best gene in population will be assigned to be elite.
   *
   */
  void __impl_select() {
    const double prevEliteFitness = _eliteIt->_Fitness;
    std::vector<GeneIt_t> iterators;
    iterators.clear();
    iterators.reserve(this->_population.size());
    auto GeneItCmp = [](GeneIt_t a, GeneIt_t b) { return isBetter(a->_Fitness, b->_Fitness); };

    for (auto it = this->_population.begin(); it != this->_population.end(); ++it) {
      iterators.emplace_back(it);
    }

    std::sort(iterators.begin(), iterators.end(), GeneItCmp);

    while (this->_population.size() > this->_option.populationSize) {
      this->_population.erase(iterators.back());
      iterators.pop_back();
    }

    GeneIt_t curBest = iterators.front();
    if (!isBetter(curBest->_Fitness, prevEliteFitness)) {
      this->_failTimes++;
      _eliteIt = curBest;
    } else {
      this->_failTimes = 0;
      _eliteIt = curBest;
    }

    this->_population.emplace_back(*_eliteIt);
  }
};

}  //    namespace heu

#endif  // HEU_SOGA_H
