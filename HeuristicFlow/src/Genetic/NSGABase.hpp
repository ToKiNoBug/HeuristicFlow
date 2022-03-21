// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_NSGABASE_HPP
#define EIGEN_HEU_NSGABASE_HPP

#include "InternalHeaderCheck.h"
#include "MOGABase.hpp"

namespace Eigen {

namespace internal {

template <typename Var_t, int ObjNum, FitnessOption fOpt, RecordOption rOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_>
class NSGABase : public MOGABase<Var_t, ObjNum, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
 private:
  using Base_t = MOGABase<Var_t, ObjNum, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  NSGABase(){};
  virtual ~NSGABase(){};

  EIGEN_HEU_MAKE_GABASE_TYPES(Base_t)
  using Fitness_t = typename Base_t::Fitness_t;

  /**
   * @brief basical unit for NS
   *
   */
  struct infoUnitBase {
   public:
    Fitness_t fitnessCache;
    /** @brief Genes in population that strong domain this gene
     */
    size_t domainedByNum;

    /** @brief iterator to related gene
     */
    GeneIt_t iterator;
  };  //  infoUnitBase

  void initializePop() {
    sortSpace.clear();
    sortSpace.reserve(2 * this->_option.populationSize);
    Base_t::initializePop();
  }

 protected:
  std::vector<infoUnitBase*> sortSpace;
  std::list<std::vector<infoUnitBase*>> pfLayers;

  static bool sortByDominatedNum(const infoUnitBase* A, const infoUnitBase* B) {
    return (A->domainedByNum) < (B->domainedByNum);
  }

  // calculate domainedByNum
  virtual void calculateDominatedNum() {
    const size_t popSizeBefore = sortSpace.size();
#ifdef EIGEN_HAS_OPENMP
    static const int32_t thN = Eigen::nbThreads();
#pragma omp parallel for schedule(dynamic, popSizeBefore / thN)
    for (int ed = 0; ed < popSizeBefore; ed++) {
      sortSpace[ed]->domainedByNum = 0;
      for (size_t er = 0; er < popSizeBefore; er++) {
        if (er == ed) continue;
        sortSpace[ed]->domainedByNum +=
            Pareto<ObjNum, fOpt>::isStrongDominate(&(sortSpace[er]->fitnessCache), &(sortSpace[ed]->fitnessCache));
      }
    }

#else
    for (size_t ed = 0; ed < popSizeBefore; ed++) {
      sortSpace[ed]->domainedByNum = 0;
      for (size_t er = 0; er < popSizeBefore; er++) {
        if (er == ed) continue;
        sortSpace[ed]->domainedByNum +=
            Pareto<ObjNum, fOpt>::isStrongDominate(&(sortSpace[er]->fitnessCache), &(sortSpace[ed]->fitnessCache));
      }
    }
#endif

  }  // calculateDominatedNum()

  virtual void divideLayers() {
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

  void updatePF(const infoUnitBase** pfs, const size_t curFrontSize) {
    this->_pfGenes.clear();
    for (size_t i = 0; i < curFrontSize; i++) {
      this->_pfGenes.emplace(&*(pfs[i]->iterator));
    }
    if (this->prevFrontSize != curFrontSize) {
      this->_failTimes = 0;
      this->prevFrontSize = curFrontSize;
    } else {
      size_t checkSum = this->makePFCheckSum();

      if (this->prevPFCheckSum == checkSum) {
        this->_failTimes++;
      } else {
        this->_failTimes = 0;
        this->prevPFCheckSum = checkSum;
      }
    }
  }  // updatePF()

};  // NSGABase

#define EIGEN_HEU_MAKE_NSGABASE_TYPES(Base_t)           \
  EIGEN_HEU_MAKE_GABASE_TYPES(Base_t)                   \
  using infoUnitBase_t = typename Base_t::infoUnitBase; \
  using Fitness_t = typename Base_t::Fitness_t;

}  //  namespace internal

}  //  namespace Eigen

#endif  //  EIGEN_HEU_NSGABASE_HPP
