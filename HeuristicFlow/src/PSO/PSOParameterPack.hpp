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

#ifndef HEU_PSOPARAMETERPACK_HPP
#define HEU_PSOPARAMETERPACK_HPP

#include <HeuristicFlow/Global>
#include "InternalHeaderCheck.h"

namespace heu {

namespace internal {

HEU_MAKE_FUNAREA(iFun, PSO)
HEU_MAKE_FUNAREA(fFun, PSO)

/**
 * \ingroup HEU_PSO
 * \class PSOParameterPack
 * \brief This class maintains a member `Args_t _args` if `Args_t` is void. Besides it defines the
 * format of initialization function and fitness function.
 *
 * \tparam Var_t Type of decision variable.
 * \tparam Fitness_t Type of fitness
 * \tparam Arg_t Type of other pesudo global parameters.
 *
 * \sa GAAbstract for its counterpart in Genetic module
 * \sa PSOParameterPack<Var_t, Fitness_t, void> Its specialization
 */
template <class Var_t, class Fitness_t, class Arg_t = void>
class PSOParameterPack {
 public:
  ~PSOParameterPack() = default;
  /// A alias of `Args_t`
  using Args_t = Arg_t;

  /**
   * \brief Type of initialization function.
   *
   * The initialization proceudre must set the position and velocity for a particle.
   *
   * This struct has a specialization for `Args_t` is `void`.
   *
   * \param pos The position to be initialized
   * \param velocity The velocity to be initialized
   * \param pMin Minimum position
   * \param pMax Maximum position
   * \param vMax Maximum velocity (Maximum absolute value of velocity)
   */
  using iFun_t = void (*)(Var_t *pos, Var_t *velocity, const Args_t *);

  /**
   * \brief Type of fitness function
   *
   * This function must compute the fitness value with the position of a particle.
   *
   * \param pos Positon value
   * \param f Fitness value to be computed.
   */
  using fFun_t = void (*)(const Var_t *pos, const Args_t *, Fitness_t *f);

  template <iFun_t _i>
  using iFunBody = typename iFunArea_PSO<Var_t *, Var_t *, const Args_t *>::template funBody<_i>;

  template <fFun_t _i>
  using fFunBody =
      typename fFunArea_PSO<const Var_t *, const Args_t *, Fitness_t *>::template funBody<_i>;

  static inline void defaultInitializeFunctionThatShouldNotBeCalled(Var_t *, Var_t *,
                                                                    const Args_t *) noexcept {
    constexpr bool
        You_think_you_called_the_default_initialize_function_of_PSO_but_it_is_reserved_as_a_symbol =
            false;
    assert(
        You_think_you_called_the_default_initialize_function_of_PSO_but_it_is_reserved_as_a_symbol);
  }

  /**
   * \brief Set the Args object
   *
   * \param a Value of args.
   */
  inline void setArgs(const Arg_t &a) noexcept { _arg = a; }

  /**
   * \brief Get the args object.
   *
   * \return const Arg_t& A const-reference to args.
   */
  inline const Arg_t &args() const noexcept { return _arg; }

  /**
   * \brief This static member denotes that `Args_t` is not `void`.
   *
   */
  static constexpr bool HasParameters = true;

 protected:
  /**
   * \brief Pseudo-global variables inside the PSO solver.
   *
   * \note This member only exists when `Arg_t` is not `void`.
   */
  Arg_t _arg;
};

/**
 * \ingroup HEU_PSO
 * \class PSOParameterPack<Var_t, Fitness_t, void>
 * \brief Specilized for PSO without args
 *
 * Most functions are same as PSOParameterPack except there is no member named `_arg`.\n Besides,
 * iFun and fFun don't have a input of type `const Args_t*`
 *
 * \sa PSOParameterPack
 */
template <class Var_t, class Fitness_t>
class PSOParameterPack<Var_t, Fitness_t, void> {
 public:
  ~PSOParameterPack() = default;
  using Args_t = void;

  using iFun_t = void (*)(Var_t *pos, Var_t *velocity);
  using fFun_t = void (*)(const Var_t *, Fitness_t *);

  template <iFun_t _i>
  using iFunBody = typename iFunArea_PSO<Var_t *, Var_t *>::template funBody<_i>;

  template <fFun_t _i>
  using fFunBody = typename fFunArea_PSO<const Var_t *, Fitness_t *>::template funBody<_i>;

  static inline void defaultInitializeFunctionThatShouldNotBeCalled(Var_t *, Var_t *) noexcept {
    constexpr bool
        You_think_you_called_the_default_initialize_function_of_PSO_but_it_is_reserved_as_a_symbol =
            false;
    assert(
        You_think_you_called_the_default_initialize_function_of_PSO_but_it_is_reserved_as_a_symbol);
  }

  static constexpr bool HasParameters = false;
};

#define HEU_MAKE_PSOPARAMETERPACK_TYPES(Base_t) \
  using Args_t = typename Base_t::Args_t;       \
  using iFun_t = typename Base_t::iFun_t;       \
  using fFun_t = typename Base_t::fFun_t;

}  //  namespace internal

}  // namespace heu

#endif  //  HEU_PSOPARAMETERPACK_HPP
