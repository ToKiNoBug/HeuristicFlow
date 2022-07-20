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
#ifndef HEU_PSO4STD_HPP
#define HEU_PSO4STD_HPP

#include <array>
#include <vector>
#include <tuple>
#include <type_traits>

#include "InternalHeaderCheck.h"
#include "PSOOption.hpp"
#include "PSOAbstrcat.hpp"

namespace heu {

namespace internal {

template <typename Var_t, FitnessOption FitnessOpt, RecordOption RecordOpt, class Arg_t,
          BoxShape BS, typename internal::PSOParameterPack<Var_t, double, Arg_t>::iFun_t _iFun_,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::fFun_t _fFun_>
class PSO4std : public internal::PSOAbstract<Var_t, double, RecordOpt, Arg_t, BS, _iFun_, _fFun_> {
  using Base_t = internal::PSOAbstract<Var_t, double, RecordOpt, Arg_t, BS, _iFun_, _fFun_>;

 public:
  PSO4std() = default;
  ~PSO4std() = default;
  HEU_MAKE_PSOABSTRACT_TYPES(Base_t)
  friend class internal::PSOAbstract<Var_t, double, DONT_RECORD_FITNESS, Arg_t, BS, _iFun_, _fFun_>;

 protected:
  /**
   * \brief To determine whether fitness a is better than b
   *
   * \param a The first input
   * \param b The second input
   * \return true A is better than B
   * \return false A is not better than B
   */
  inline static bool isBetterThan(double a, double b) noexcept {
    if (FitnessOpt == FitnessOption::FITNESS_GREATER_BETTER) {
      return a > b;
    } else {
      return a < b;
    }
  }

  /**
   * \brief Update gBest and pBest
   *
   */
  void __impl_updatePGBest() noexcept {
    // gBest for current generation
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
   * \brief Update the position and velocity of all particles
   *
   */
  void __impl_updatePopulation() noexcept {
#ifdef HEU_HAS_OPENMP
    static const int32_t thN = threadNum();
#pragma omp parallel for schedule(dynamic, this->_population.size() / thN)
#endif  //  HEU_HAS_OPENMP
    for (int index = 0; index < this->_population.size(); index++) {
      Particle_t& i = this->_population[index];
      const double rndP = randD();
      const double rndG = randD();
      for (int idx = 0; idx < this->dimensions(); idx++) {
        i.velocity[idx] =
            this->_option.inertiaFactor * i.velocity[idx] +
            this->_option.learnFactorP * rndP * (i.pBest.position[idx] - i.position[idx]) +
            this->_option.learnFactorG * rndG * (this->gBest.position[idx] - i.position[idx]);
        if (std::abs(i.velocity[idx]) > this->_velocityMax[idx]) {
          i.velocity[idx] = sign(i.velocity[idx]) * this->_velocityMax[idx];
        }
        i.position[idx] += i.velocity[idx];
        i.position[idx] = std::max(i.position[idx], this->_posMin[idx]);
        i.position[idx] = std::min(i.position[idx], this->_posMax[idx]);
      }
    }
  }

 private:
  static_assert(!(std::is_scalar<Var_t>::value), "Var_t should be a non-scalar type");
  // static_assert(isEigenTypes == false, "Wrong specialization of PSO");
};

}  //  namespace internal

}  //  namespace heu

#endif  //  HEU_PSO4STD_HPP
