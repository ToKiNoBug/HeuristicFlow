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

#ifndef HEU_NSGABASE_HPP
#define HEU_NSGABASE_HPP

#include "InternalHeaderCheck.h"
#include "MOGABase.hpp"

namespace heu {

namespace internal {

/**
 * \ingroup HEU_GENETIC
 * \class NSGABase
 * \brief Base class for nondominated sorting series of genetic algorithms, like NSGA2 and NSGA3.
 *
 * This class implements nondominated sorting, which is common in NSGA2 and NSGA3.
 *
 * \tparam Var_t
 * \tparam ObjNum
 * \tparam fOpt
 * \tparam rOpt
 * \tparam Args_t
 * \tparam _iFun_
 * \tparam _fFun_
 * \tparam _cFun_
 * \tparam _mFun_
 */
template <typename Var_t, int ObjNum, FitnessOption fOpt, RecordOption rOpt, class Gene,
          class Args_t,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_>
class NSGABase
    : public MOGABase<Var_t, ObjNum, fOpt, rOpt, Gene, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
 private:
  using Base_t = MOGABase<Var_t, ObjNum, fOpt, rOpt, Gene, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  ~NSGABase() = default;
  HEU_MAKE_GABASE_TYPES(Base_t)
  using Fitness_t = typename Base_t::Fitness_t;

  /**
   * \brief Basical unit for NS
   *
   * This struct stores basic informations of a gene used in nondominated sorting.
   */
  struct infoUnitBase {
   public:
    /**
     * \brief A copy of fitness to boost nondominated sorting.
     *
     */
    Fitness_t fitnessCache;

    /**
     * \brief Genes in population that strong dominate this gene
     *
     * Value 0 means it's a member of the pareto
     * front.
     */
    size_t domainedByNum;

    /**
     * \brief Iterator to related gene
     */
    GeneIt_t iterator;

  };  //  infoUnitBase

  /**
   * \brief Reimplemented
   *
   */
  inline void initializePop() noexcept {
    sortSpace.clear();
    sortSpace.reserve(2 * this->_option.populationSize);
    Base_t::initializePop();
  }

 protected:
  /**
   * \brief A vector of infoUnitBase*. This vector is where NS-sorting will take place since copying
   * swaping pointers is much faster than swaping structs.
   *
   */
  std::vector<infoUnitBase*> sortSpace;

  /**
   * \brief A list of vectors that divides the current population into several nondominated layers.
   *
   */
  std::list<std::vector<infoUnitBase*>> pfLayers;

  /**
   * \brief Compare function used to sort by dominatedNum in ascend.
   *
   * \param A The first input
   * \param B The second input
   * \return true A is dominated by less genes than B
   * \return false A is dominated by more or equal genes than B
   */
  inline static bool sortByDominatedNum(const infoUnitBase* A, const infoUnitBase* B) noexcept {
    return (A->domainedByNum) < (B->domainedByNum);
  }

  // calculate domainedByNum

  /**
   * \brief Function to calculate infoUnitBase::domainedByNum of a population.
   *
   */
  void calculateDominatedNum() noexcept {
    const size_t popSizeBefore = sortSpace.size();
#ifdef HEU_HAS_OPENMP
    static const int32_t thN = threadNum();
#pragma omp parallel for schedule(dynamic, popSizeBefore / thN)
    for (int ed = 0; ed < popSizeBefore; ed++) {
      sortSpace[ed]->domainedByNum = 0;
      for (size_t er = 0; er < popSizeBefore; er++) {
        if (er == ed) continue;
        sortSpace[ed]->domainedByNum += Pareto<ObjNum, fOpt>::isStrongDominate(
            &(sortSpace[er]->fitnessCache), &(sortSpace[ed]->fitnessCache));
      }
    }

#else
    for (size_t ed = 0; ed < popSizeBefore; ed++) {
      sortSpace[ed]->domainedByNum = 0;
      for (size_t er = 0; er < popSizeBefore; er++) {
        if (er == ed) continue;
        sortSpace[ed]->domainedByNum += Pareto<ObjNum, fOpt>::isStrongDominate(
            &(sortSpace[er]->fitnessCache), &(sortSpace[ed]->fitnessCache));
      }
    }
#endif

  }  // calculateDominatedNum

  /**
   * \brief Divide the sorted population (is sortSpace) into a few non-dominated layers.
   *
   */
  void divideLayers() noexcept {
    std::sort(sortSpace.begin(), sortSpace.end(), sortByDominatedNum);
    pfLayers.clear();
    const size_t popSizeBef = sortSpace.size();
    size_t curDM = -1;
    for (auto i : sortSpace) {
      if (curDM != i->domainedByNum) {
        curDM = i->domainedByNum;
        pfLayers.emplace_back();
        pfLayers.back().reserve(popSizeBef);
      }
      pfLayers.back().emplace_back(i);
    }
  }

  /**
   * \brief Update the pareto front
   *
   */
  void updatePF(const infoUnitBase** pfs, const size_t curFrontSize) noexcept {
    this->_pfGenes.clear();
    for (size_t i = 0; i < curFrontSize; i++) {
      this->_pfGenes.emplace(&*(pfs[i]->iterator));
    }
  }  // updatePF()

};  // NSGABase

#define HEU_MAKE_NSGABASE_TYPES(Base_t)                 \
  HEU_MAKE_GABASE_TYPES(Base_t)                         \
  using infoUnitBase_t = typename Base_t::infoUnitBase; \
  using Fitness_t = typename Base_t::Fitness_t;

}  //  namespace internal

}  //  namespace heu

#endif  //  HEU_NSGABASE_HPP
