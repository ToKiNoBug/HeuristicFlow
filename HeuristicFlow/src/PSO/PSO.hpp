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
#include "PSOAbstrcat.hpp"

#include "PSO4Eigen.hpp"
#include "PSO4std.hpp"

namespace heu {

namespace internal {}  //  namespace internal

/**
 * \ingroup HEU_PSO
 * \class PSO
 * \brief Generalized PSO solver.
 *
 * All default value for template parameters are listed in braces
 *
 * \tparam Var_t Type of decision variable.
 * \tparam BS The shape of box-constraint of your PSO solver.
 * \tparam FitnessOpt Trainning direction (FITNESS_LESS_BETTER)
 * \tparam RecordOpt Record trainning curve or not. (DONT_RECORD_FITNESS)
 * \tparam Arg_t Pseudo-global other args stored in the solver. (void)
 * \tparam _fFun_ Fitness function at compile time (nullptr)
 * \tparam _iFun_ Initialization function at compile time (A internal default initializer).
 *
 * \sa PSOOption
 * \sa SOGA for the meanning of `Arg_t`
 * \sa BoxShape
 *
 * PSO solvers have different APIs for different template parameters.
 *
 * ## These APIs are common:
 * - `void setOption(const PSOOption&)` to set the option of a PSO solver.
 * - `void initializePop()` to initialize the whole population.
 * - `void run()` to run the PSO algorithm.
 * - `double bestFitness() const` returns the best fitness value.
 * - `const PSOOption& option() const` returns a const-ref to the PSOOption object.
 * - `size_t generation() const` returns the generations that the solver has spend.
 * - `size_t failTimes() const` returns the failtimes for current population.
 * - `const auto & posMin() const` returns a const-ref to the minimum value of positoin.
 * - `const auto & posMax() const` returns a const-ref to the maximum value of position.
 * - `const auto& velocityMax() const` returns a const-ref to the maximum value of the abstract
 * value of velocity.
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
 * ## These APIs exist when BS is SQUARE_BOX :
 * - `void setPVRange(const Var_t& pMin, const Var_t& pMax, const Var_t& vMax)` to set the value of
 * posMin, posMax, and velocityMax.
 * - `void setRange(const Var_t & pMin, const Var_t & pMax)`
 * - `void setMaxVelocity(const Var_t & maxVelocity)`
 *
 * ### Otherwise it will be :
 * - `void setPVRange(const Scalar_t pMin, const Scalar_t pMax,const Scalar_t vMax)` to set the
 * value of posMin, posMax, and velocityMax.
 * - `void setRange(const Scalar_t  pMin, const Scalar_t  pMax)`
 * - `void setMaxVelocity(const Scalar_t  maxVelocity)`
 *
 * ## These APIs exist when `Args_t` is not `void` :
 * - `const Arg_t &args() const` returns a const-ref to args in the solver.
 * - `void setArgs(const Arg_t &)` sets the value of args.
 *
 *
 * ## These APIs exist when template parameter `_iFun_` is `nullptr` :
 * - `setiFun(iFun_t)` sets the initialization function.
 *
 * ## These APIs exist when template parameter `_fFun_` is `nullptr` :
 * - `setfFun(fFun_t)` sets the fitness function.
 *
 *
 * ## These APIs exist when `Var_t` CAN be a matrix :
 * - `int boxRows() const` return the rows of posMin, posMax and velocityMax.
 * - `int boxCols() const` return the cols of posMin, posMax and velocityMax.
 *
 * ## These APIs exist when `Var_t` has a dynamic size AND CAN be a matrix :
 * - `void setDimensions(int r,int c)` sets the size of posMin, posMax and velocityMax to (r,c).
 *
 * ## These APIs exist when `Var_t` has a dynamic size AND MUST be a vector :
 * - `void setDimensions(int d)` sets the size of posMin, posMax and velocityMax to (d,1).
 *
 *
 *
 */
template <typename Var_t, BoxShape BS = BoxShape::RECTANGLE_BOX,
          FitnessOption FitnessOpt = FITNESS_LESS_BETTER,
          RecordOption RecordOpt = DONT_RECORD_FITNESS, class Arg_t = void,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::fFun_t _fFun_ = nullptr,
          typename internal::PSOParameterPack<Var_t, double, Arg_t>::iFun_t _iFun_ =
              internal::PSOParameterPack<Var_t, double,
                                         Arg_t>::defaultInitializeFunctionThatShouldNotBeCalled>
class PSO
    : public std::conditional<
          isEigenClass<Var_t>::value,
          typename internal::PSO4Eigen<Var_t, FitnessOpt, RecordOpt, Arg_t, BS, _iFun_, _fFun_>,
          typename internal::PSO4std<Var_t, FitnessOpt, RecordOpt, Arg_t, BS, _iFun_,
                                     _fFun_>>::type {
  using Base_t = typename std::conditional<
      isEigenClass<Var_t>::value,
      typename internal::PSO4Eigen<Var_t, FitnessOpt, RecordOpt, Arg_t, BS, _iFun_, _fFun_>,
      typename internal::PSO4std<Var_t, FitnessOpt, RecordOpt, Arg_t, BS, _iFun_, _fFun_>>::type;

 public:
  PSO() = default;
  ~PSO() = default;

  HEU_RELOAD_MEMBERFUCTION_RUN
  HEU_MAKE_PSOABSTRACT_TYPES(Base_t)

  /**
   * \brief Function used to provide a result for recording
   *
   * \return double The best fitness to be recorded
   */
  inline double bestFitness() const noexcept { return this->gBest.fitness; }
};
}  // namespace heu

#endif  // HEU_PSO_HPP
