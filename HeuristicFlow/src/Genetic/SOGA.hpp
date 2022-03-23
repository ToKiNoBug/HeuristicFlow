// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_SOGA_H
#define EIGEN_HEU_SOGA_H

#include "InternalHeaderCheck.h"
#include "GABase.hpp"

namespace Eigen {

/**
 * \brief Single-object genetic solver.
 *
 *  @tparam Var_t  Type of decisition variable.
 *  @tparam fOpt Whether greater fitness value means better. SOGA will always try to find the best, so it's vital to
 * tell SOGA which direction is right.
 *  @tparam Record  Whether the solver records fitness changelog.
 *  @tparam Args_t  Type of other parameters.
 * \tparam _iFun_ Function to initialize an individual in population. This function will be called only when
 * initializing. Use nullptr if you hope to determine it at runtime. nullptr as default value.
 * \tparam _fFun_ Funtion to compute fitness for any
 * individual. Use nullptr if you hope to determine it at runtime. nullptr as default value.
 * \tparam _cFun_ Function to apply crossover. Use
 * nullptr if you hope to determine it at runtime. nullptr as default value.
 * \tparam _mFun_ Function to apply mutation. Use nullptr if you hope to
 * determine it at runtime. nullptr as default value.
 */
template <typename Var_t, FitnessOption fOpt = FITNESS_LESS_BETTER, RecordOption Record = DONT_RECORD_FITNESS,
          class Args_t = void, typename internal::GAAbstract<Var_t, double, Args_t>::initializeFun _iFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::fitnessFun _fFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::crossoverFun _cFun_ = nullptr,
          typename internal::GAAbstract<Var_t, double, Args_t>::mutateFun _mFun_ = nullptr>
class SOGA : public internal::GABase<Var_t, double, Record, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
 private:
  using Base_t = internal::GABase<Var_t, double, Record, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  EIGEN_HEU_MAKE_GABASE_TYPES(Base_t)

  /**
   * \brief Return the best fitness of current population
   *
   * \return double best fitness value.
   */
  virtual double bestFitness() const { return _eliteIt->_Fitness; }

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
  virtual void select() {
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

}  //    namespace Eigen

#endif  // EIGEN_HEU_SOGA_H
