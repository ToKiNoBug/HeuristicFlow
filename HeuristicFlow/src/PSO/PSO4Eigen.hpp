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
#ifndef HEU_PSO4EIGEN_HPP
#define HEU_PSO4EIGEN_HPP

#include <array>
#include <vector>
#include <tuple>
#include <type_traits>

#include "InternalHeaderCheck.h"
#include "PSOOption.hpp"
#include "PSOBase.hpp"

namespace heu {

namespace internal {
// Partial specilization for PSO using Eigen's fix-sized Array
template <typename Var_t, int DIM, FitnessOption FitnessOpt, RecordOption RecordOpt, class Arg_t,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::iFun_t _iFun_,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::fFun_t _fFun_>
class PSO4Eigen : public internal::PSOBase<Var_t, DIM, double, RecordOpt, Arg_t, _iFun_, _fFun_> {
  using Base_t = internal::PSOBase<Var_t, DIM, double, RecordOpt, Arg_t, _iFun_, _fFun_>;

 public:
  HEU_MAKE_PSOABSTRACT_TYPES(Base_t)
  friend class internal::PSOAbstract<Var_t, double, DONT_RECORD_FITNESS, Arg_t, _iFun_, _fFun_>;

 protected:
  static bool isBetterThan(double a, double b) {
    if (FitnessOpt == FitnessOption::FITNESS_GREATER_BETTER) {
      return a > b;
    } else {
      return a < b;
    }
  }

  /**
   * \brief Update the value of pBest and gBest
   *
   */
  void __impl_updatePGBest() {
    Point_t* curGBest = &this->_population.front().pBest;

    for (Particle_t& i : this->_population) {
      if (isBetterThan(i.fitness, i.pBest.fitness)) {
        i.pBest = i;
      }

      if (isBetterThan(i.pBest.fitness, curGBest->fitness)) {
        curGBest = &i.pBest;
      }
    }

    if (isBetterThan(curGBest->fitness, this->gBest.fitness)) {
      this->_failTimes = 0;
      this->gBest = *curGBest;
    } else {
      this->_failTimes++;
    }
  }

  /**
   * \brief Update the position and velocity of each particle
   *
   */
  void __impl_updatePopulation() {
#ifdef HEU_HAS_OPENMP
    static const int32_t thN = threadNum();
#pragma omp parallel for schedule(dynamic, this->_population.size() / thN)
#endif  //  HEU_HAS_OPENMP
    for (int idx = 0; idx < (int)this->_population.size(); idx++) {
      Particle_t& i = this->_population[idx];
      const Scalar_t lFP = randD(), lFG = randD();
      i.velocity = this->_option.inertiaFactor * i.velocity +
                   this->_option.learnFactorP * lFP * (i.pBest.position - i.position) +
                   this->_option.learnFactorG * lFG * (this->gBest.position - i.position);

      i.velocity = i.velocity.min(this->_velocityMax);
      i.velocity = i.velocity.max(-this->_velocityMax);

      i.position += i.velocity;

      i.position = i.position.min(this->_posMax);
      i.position = i.position.max(this->_posMin);
    }
  }

 private:
  static_assert(DIM > 0 || DIM == Eigen::Dynamic, "Invalid template parameter DIM");
};

}  //  namespace internal

}  //  namespace heu

#endif  //  HEU_PSO4EIGEN_HPP
