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

#ifndef HEU_MOGAABSTRACT_HPP
#define HEU_MOGAABSTRACT_HPP

#include <queue>
#include <functional>
#include <unordered_set>

#include <HeuristicFlow/Global>
#include <HeuristicFlow/EAGlobal>
#include "InternalHeaderCheck.h"
#include "GABase.hpp"

namespace heu {

namespace internal {

/**
 * \ingroup HEU_GENETIC
 * \class MOGAAbstract
 * \brief Base class for multi-objective genetic algorithm solver.
 *
 * This class maintains a pareto frontier.
 *
 * \tparam Var_t
 * \tparam ObjNum Number of objectives. Use Eigen::Dynamic for runtime determined objective numbers.
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
class MOGAAbstract : public GABase<Var_t, Eigen::Array<double, ObjNum, 1>, rOpt, Gene, Args_t,
                                   _iFun_, _fFun_, _cFun_, _mFun_> {
 private:
  using Base_t = GABase<Var_t, Eigen::Array<double, ObjNum, 1>, rOpt, Gene, Args_t, _iFun_, _fFun_,
                        _cFun_, _mFun_>;
  static_assert(ObjNum > 0 || ObjNum == Eigen::Dynamic, "Invalid template parameter Dim");
  static_assert(ObjNum != 1, "You assigend 1 objective in multi-objective problems");

 public:
  ~MOGAAbstract() = default;
  HEU_MAKE_GABASE_TYPES(Base_t)

  /**
   * \brief Type of fitness is stored in an Eigen vector of double.
   *
   */
  using Fitness_t = Eigen::Array<double, ObjNum, 1>;

  /// get pareto front in vec

  /**
   * \brief Get Pareto front consisted of a vector of Fitness_t.
   *
   * \param front vector of fitness.
   */
  inline void paretoFront(std::vector<Fitness_t>& front) const noexcept {
    front.clear();
    front.reserve(_pfGenes.size());
    for (const Gene* i : _pfGenes) {
      front.emplace_back(i->_Fitness);
    }
    return;
  }

  /**
   * \brief Get Pareto front consists of a vector of std::pair<const Var_t*,const Fitness_t*>. It
   * provides both solution(Var_t) and objective value(Fitness_t).
   *
   * \param front vector of Var and fitness.
   */
  inline void paretoFront(
      std::vector<std::pair<const Var_t*, const Fitness_t*>>& front) const noexcept {
    front.clear();
    front.reserve(_pfGenes.size());
    for (const Gene* i : _pfGenes) {
      front.emplace_back(std::make_pair(&(i->self), &(i->_Fitness)));
    }
  }

  /**
   * \brief Returns a const reference to member _pfGenes.
   *
   * \return const std::unordered_set<const Gene*>& A const reference to PF.
   */
  inline const std::unordered_set<const Gene*>& pfGenes() const noexcept { return _pfGenes; }

  /**
   * \brief This function is reimplemented to resize the PF.
   *
   */
  inline void initializePop() noexcept {
    this->_pfGenes.clear();
    this->_pfGenes.reserve(this->_option.populationSize * 2);
    Base_t::initializePop();
  }

  /**
   * \brief Compute ideal point
   *
   * An ideal point is a fitness whose componment of each objecvite is the minimium value of the
   * whole population.
   *
   * \note The corresponding decision variable of the ideal point usually doesn't exist.
   *
   * \return Fitness_t ideal point
   */
  Fitness_t bestFitness() const noexcept {
    Fitness_t best = this->_population.front()._Fitness;
    for (const Gene& i : this->_population) {
      if (fOpt == FitnessOption::FITNESS_GREATER_BETTER) {
        best = best.max(i._Fitness);
      } else {
        best = best.min(i._Fitness);
      }
    }
    return best;
  }

 protected:
  std::unordered_set<const Gene*> _pfGenes;  ///< A hash set to store the whole PF

  /*
   * \brief Compute the hash checksum of current PF
   *
   * PF checksum is computed with the hash of addresses of every Gene that consist of PF.
   *
   * Genes are sorted by their address, and each hash of address are sumed up by xor.
   *
   * This method works because elitisim strategy is applied. Actually the value of a gene won't
   * change once it's initialized(created during initialization, crossover or mutation). And a gene
   * will never be simply copied. Thus hashing the address is fast and it works.
   *
   * \return size_t The checksum of PF
   */
  /*
  size_t makePFCheckSum() const noexcept {
    std::vector<const Gene*> pfvec;
    pfvec.reserve(_pfGenes.size());
    for (auto i : _pfGenes) {
      pfvec.emplace_back(i);
    }
    std::sort(pfvec.begin(), pfvec.end());

    size_t checkSum = std::hash<const void*>()(pfvec[0]);
    for (size_t i = 1; i < pfvec.size(); i++) {
      checkSum ^= std::hash<const void*>()(pfvec[i]);
    }
    return checkSum;
  }*/

};  // MOGAAbstract

}  //  namespace internal

}  //  namespace heu

#endif  //   HEU_MOGAABSTRACT_HPP
