// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_PSOPARAMETERPACK_HPP
#define EIGEN_HEU_PSOPARAMETERPACK_HPP

#include "../../Global"
#include "InternalHeaderCheck.h"

namespace Eigen {

namespace internal {

EIGEN_HEU_MAKE_FUNAREA(iFun, iFun, PSO)
EIGEN_HEU_MAKE_FUNAREA(fFun, fFun, PSO)

/**
 * \ingroup HEU_PSO
 * \class PSOParameterPack
 * \brief This class maintains a member `Args_t _args` if `Args_t` is void. Besides it defines the format of
 * initialization function and fitness function.
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
  using iFun_t = void (*)(Var_t *pos, Var_t *velocity, const Var_t *pMin, const Var_t *pMax, const Var_t *vMax,
                          const Args_t *);

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
  using iFunBody = typename iFunArea_PSO<void, Var_t *, Var_t *, const Var_t *, const Var_t *, const Var_t *,
                                         const Args_t *>::template funBody<_i>;

  template <fFun_t _i>
  using fFunBody = typename fFunArea_PSO<void, const Var_t *, const Args_t *, Fitness_t *>::template funBody<_i>;

  /**
   * \brief Set the Args object
   *
   * \param a Value of args.
   */
  void setArgs(const Arg_t &a) { _arg = a; }

  /**
   * \brief Get the args object.
   *
   * \return const Arg_t& A const-reference to args.
   */
  const Arg_t &args() const { return _arg; }

  /**
   * \brief This static member denotes that `Args_t` is not `void`.
   *
   */
  static const bool HasParameters = true;

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
 * \brief
 *
 */
template <class Var_t, class Fitness_t>
class PSOParameterPack<Var_t, Fitness_t, void> {
 public:
  using Args_t = void;

  using iFun_t = void (*)(Var_t *pos, Var_t *velocity, const Var_t *pMin, const Var_t *pMax, const Var_t *vMax);
  using fFun_t = void (*)(const Var_t *, Fitness_t *);

  template <iFun_t _i>
  using iFunBody =
      typename iFunArea_PSO<void, Var_t *, Var_t *, const Var_t *, const Var_t *, const Var_t *>::template funBody<_i>;

  template <fFun_t _i>
  using fFunBody = typename fFunArea_PSO<void, const Var_t *, Fitness_t *>::template funBody<_i>;

  static const bool HasParameters = false;
};

#define EIGEN_HEU_MAKE_PSOPARAMETERPACK_TYPES(Base_t) \
  using Args_t = typename Base_t::Args_t;             \
  using iFun_t = typename Base_t::iFun_t;             \
  using fFun_t = typename Base_t::fFun_t;

}  //  namespace internal

}  // namespace Eigen

#endif  //  EIGEN_HEU_PSOPARAMETERPACK_HPP
