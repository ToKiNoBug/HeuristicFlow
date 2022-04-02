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

#ifndef HEU_PSO_HPP
#define HEU_PSO_HPP
#include <array>
#include <vector>
#include <tuple>
#include <type_traits>

#include "InternalHeaderCheck.h"
#include "PSOOption.hpp"
#include "PSOBase.hpp"

namespace heu {

/**
 * \ingroup HEU_PSO
 * \class PSO
 * \brief Generalized PSO solver.
 *
 * All default value for template parameters are listed in braces
 *
 * \tparam Var_t Type of decision variable.
 * \tparam DIM Dimensional of Var_t
 * \tparam isEigenTypes If Var_t is a Eigen's Array/Matrix (true)
 * \tparam FitnessOpt Trainning direction (FITNESS_LESS_BETTER)
 * \tparam RecordOpt Record trainning curve or not. (DONT_RECORD_FITNESS)
 * \tparam Arg_t Pseudo-global other args stored in the solver. (void)
 * \tparam _iFun_ Initialization function at compile time (nullptr)
 * \tparam _fFun_ Fitness function at compile time (nullptr)
 *
 * \sa SOGA for the meanning of `Arg_t`
 *
 * \note This class implements PSO for the condition that `isEigenTypes` is `false`.
 *
 * PSO solvers have different APIs for different template parameters.
 *
 * ## These APIs are common:
 * - `void setOption(const PSOOption&)` to set the option of a PSO solver.
 * - `void setPVRange(const Var_t& pMin, const Var_t& pMax, const Var_t& vMax)` to set the value of posMin, posMax, and
 * velocityMax.
 * - `void setPVRange(double pMin, double pMax, double vMax)` will also set these value but it shapes the box into a
 * square box.
 * - `void initializePop()` to initialize the whole population.
 * - `void run()` to run the PSO algorithm.
 * - `double bestFitness() const` returns the best fitness value.
 * - `const PSOOption& option() const` returns a const-ref to the PSOOption object.
 * - `size_t generation() const` returns the generations that the solver has spend.
 * - `size_t failTimes() const` returns the failtimes for current population.
 * - `const Var_t& posMin() const` returns a const-ref to the minimum value of positoin.
 * - `const Var_t& posMax() const` returns a const-ref to the maximum value of position.
 * - `const Var_t& velocityMax() const` returns a const-ref to the maximum value of the abstract value of velocity.
 * - `const std::vector<Particle>& population() const` returns a const-ref to the population.
 * - `const Point& globalBest() const` returns a const-ref to the best solution that has ever found.
 * - `int dimensions() const` returns the dimensions of decision variables.
 * - `solver_t::HasParameters` marks whether `Arg_t` is `void`
 * - `typename solver_t::iFun_t` is the type of initialization function
 * - `typename solver_t::fFun_t` is the type of fitness function.
 * - `iFun_t iFun() const` returns the initialization function.
 * - `fFun_t fFun() const` returns the fitness function.
 *
 *
 * ## These APIs exist when `Args_t` is not `void` :
 * - `const Arg_t &args() const` returns a const-ref to args in the solver.
 * - `void setArgs(const Arg_t &)` set the value of args.
 *
 *
 * ## These APIs exist when template parameter `_iFun_` is `nullptr` :
 * - `setiFun(iFun_t)` sets the initialization function.
 *
 * ## These APIs exist when template parameter `_fFun_` is `nullptr` :
 * - `setfFun(fFun_t)` sets the fitness function.
 *
 * ## These APIs exist when template parameter `ObjNum` is `Eigen::Dynamic` :
 * - `void setDimensions(int d)` set the size of posMin, posMax and velocityMax to d.
 *
 */
template <typename Var_t, int DIM, bool isEigenTypes = true, FitnessOption FitnessOpt = FITNESS_LESS_BETTER,
          RecordOption RecordOpt = DONT_RECORD_FITNESS, class Arg_t = void,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::iFun_t _iFun_ = nullptr,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::fFun_t _fFun_ = nullptr>
class PSO : public internal::PSOBase<Var_t, DIM, double, RecordOpt, Arg_t, _iFun_, _fFun_> {
  using Base_t = internal::PSOBase<Var_t, DIM, double, RecordOpt, Arg_t, _iFun_, _fFun_>;

 public:
  PSO() {}
  ~PSO() {}
  HEU_MAKE_PSOABSTRACT_TYPES(Base_t)
  friend class internal::PSOAbstract<Var_t, double, DONT_RECORD_FITNESS, Arg_t, _iFun_, _fFun_>;

  static const ContainerOption Flag =
      (std::is_same<Var_t, stdVecD_t<DIM>>::value) ? ContainerOption::Std : ContainerOption::Custom;

  /**
   * \brief Get the current gBest fitness
   *
   * \return double Current global best fitness
   */
  double bestFitness() const { return this->gBest.fitness; }

  /**
   * \brief Run the PSO solver
   *
   */
  inline void run() { this->template __impl_run<PSO>(); }

 protected:
  /**
   * \brief To determine whether fitness a is better than b
   *
   * \param a The first input
   * \param b The second input
   * \return true A is better than B
   * \return false A is not better than B
   */
  inline static bool isBetterThan(double a, double b) {
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
  void __impl_updatePGBest() {
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
  void __impl_updatePopulation() {
#ifdef EIGEN_HAS_OPENMP
    static const int32_t thN = Eigen::nbThreads();
#pragma omp parallel for schedule(dynamic, this->_population.size() / thN)
#endif  //  EIGEN_HAS_OPENMP
    for (int index = 0; index < this->_population.size(); index++) {
      Particle_t& i = this->_population[index];
      const double rndP = randD();
      const double rndG = randD();
      for (int idx = 0; idx < this->dimensions(); idx++) {
        i.velocity[idx] = this->_option.inertiaFactor * i.velocity[idx] +
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
  static_assert(DIM != 0, "Template parameter DIM cannot be 0. For dynamic dims, use Eigen::Dynamic");
  static_assert(DIM > 0 || DIM == Eigen::Dynamic, "Invalid template parameter DIM");
  static_assert(isEigenTypes == false, "Wrong specialization of PSO");
};

//

/**
 * \ingroup HEU_PSO
 * \brief Convenient typedef for stdArray (fix-sized and Runtime sized)
 *
 * \tparam DIM Dimensions of decision variable. Use `Eigen::Dynamic` for runtime determined.
 * \tparam FitnessOpt Optimization direction
 * \tparam RecordOpt Whether to record the fitness of each generation when running.
 * \tparam Arg_t Pseudo-global other args stored in the solver. (void)
 * \tparam _iFun_ Initialization function at compile time (nullptr)
 * \tparam _fFun_ Fitness function at compile time (nullptr)
 */
template <int DIM, FitnessOption FitnessOpt, RecordOption RecordOpt, class Arg_t = void,
          typename internal::PSOParameterPack<stdVecD_t<DIM>, double, Arg_t>::iFun_t _iFun_ = nullptr,
          typename internal::PSOParameterPack<stdVecD_t<DIM>, double, Arg_t>::fFun_t _fFun_ = nullptr>
using PSO_std = PSO<stdVecD_t<DIM>, DIM, false, FitnessOpt, RecordOpt, Arg_t, _iFun_, _fFun_>;

/**
 * \ingroup HEU_PSO
 * \brief Convenient typedef for Eigen Arrays (fix-sized and Runtime sized)
 *
 * \tparam DIM Dimensions of decision variable. Use `Eigen::Dynamic` for runtime determined.
 * \tparam FitnessOpt Optimization direction
 * \tparam RecordOpt Whether to record the fitness of each generation when running.
 * \tparam Arg_t Pseudo-global other args stored in the solver. (void)
 * \tparam _iFun_ Initialization function at compile time (nullptr)
 * \tparam _fFun_ Fitness function at compile time (nullptr)
 */
template <int DIM, FitnessOption FitnessOpt, RecordOption RecordOpt, class Arg_t = void,
          typename internal::PSOParameterPack<Eigen::Array<double, DIM, 1>, double, Arg_t>::iFun_t _iFun_ = nullptr,
          typename internal::PSOParameterPack<Eigen::Array<double, DIM, 1>, double, Arg_t>::fFun_t _fFun_ = nullptr>
using PSO_Eigen = PSO<Eigen::Array<double, DIM, 1>, DIM, true, FitnessOpt, RecordOpt, Arg_t, _iFun_, _fFun_>;

/**
 * \ingroup HEU_PSO
 * \class PSO<Var_t, DIM, true, FitnessOpt, RecordOpt, Arg_t, _iFun_, _fFun_>
 * \brief Partial specilization for PSO using Eigen's fix-sized Array
 *
 * This specialization has the same function with PSO, but uses Eigen's api as much as possible, which provides chances
 * to be boosted.
 *
 * This class has exactly same API with PSO.
 *
 * \sa PSO For API format
 *
 * \tparam Var_t
 * \tparam DIM
 * \tparam FitnessOpt
 * \tparam RecordOpt
 * \tparam Arg_t
 * \tparam _iFun_
 * \tparam _fFun_
 */
template <typename Var_t, int DIM, FitnessOption FitnessOpt, RecordOption RecordOpt, class Arg_t,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::iFun_t _iFun_,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::fFun_t _fFun_>
class PSO<Var_t, DIM, true, FitnessOpt, RecordOpt, Arg_t, _iFun_, _fFun_>
    // Partial specilization for PSO using Eigen's fix-sized Array
    : public internal::PSOBase<Var_t, DIM, double, RecordOpt, Arg_t, _iFun_, _fFun_> {
  using Base_t = internal::PSOBase<Var_t, DIM, double, RecordOpt, Arg_t, _iFun_, _fFun_>;

 public:
  HEU_MAKE_PSOABSTRACT_TYPES(Base_t)
  friend class internal::PSOAbstract<Var_t, double, DONT_RECORD_FITNESS, Arg_t, _iFun_, _fFun_>;

  /**
   * \brief Function used to provide a result for recording
   *
   * \return double The best fitness to be recorded
   */
  double bestFitness() const { return this->gBest.fitness; }

  /**
   * \brief Run the PSO solver
   *
   */
  inline void run() { this->template __impl_run<PSO>(); }

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
#ifdef EIGEN_HAS_OPENMP
    static const int32_t thN = Eigen::nbThreads();
#pragma omp parallel for schedule(dynamic, this->_population.size() / thN)
#endif  //  EIGEN_HAS_OPENMP
    for (int idx = 0; idx < (int)this->_population.size(); idx++) {
      Particle_t& i = this->_population[idx];
      i.velocity = this->_option.inertiaFactor * i.velocity +
                   this->_option.learnFactorP * randD() * (i.pBest.position - i.position) +
                   this->_option.learnFactorG * randD() * (this->gBest.position - i.position);

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

}  // namespace heu

#endif  // HEU_PSO_HPP
