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

#ifndef HEU_AOSPARAMETERPACK_HPP
#define HEU_AOSPARAMETERPACK_HPP

#include <HeuristicFlow/Global>

#include "InternalHeaderCheck.hpp"

namespace heu {

namespace internal {

HEU_MAKE_FUNAREA(iFun, AOS)

HEU_MAKE_FUNAREA(fFun, AOS)

template <typename Var_t, class Fitness_t, class Arg_t>
class AOSParameterPack {
 public:
  using Args_t = Arg_t;
  using iFun_t = void (*)(Var_t *, const Arg_t *);
  using fFun_t = void (*)(const Var_t *, const Arg_t *, Fitness_t *);

  template <iFun_t _i>
  using iFunBody = typename iFunArea_AOS<Var_t *, const Arg_t *>::template funBody<_i>;
  template <fFun_t _f>
  using fFunBody =
      typename fFunArea_AOS<const Var_t *, const Arg_t *, Fitness_t *>::template funBody<_f>;

  static void defaultInitializeFunctionThatShouldNotBeCalled(Var_t *, const Arg_t *) {
    constexpr bool
        You_think_you_called_the_default_initialize_function_of_AOS_but_it_is_reserved_as_a_symbol =
            false;
    assert(
        You_think_you_called_the_default_initialize_function_of_AOS_but_it_is_reserved_as_a_symbol);
  }

  static constexpr bool hasParameters = true;

 protected:
  Arg_t _arg;

 public:
  inline const Arg_t &arg() const noexcept { return _arg; }
  inline Arg_t &arg() noexcept { return _arg; }
  inline void setArg(const Arg_t &_a) noexcept { _arg = _a; }
};

template <typename Var_t, class Fitness_t>
class AOSParameterPack<Var_t, Fitness_t, void> {
 public:
  using Args_t = void;
  using iFun_t = void (*)(Var_t *);
  using fFun_t = void (*)(const Var_t *, Fitness_t *);

  template <iFun_t _i>
  using iFunBody = typename iFunArea_AOS<Var_t *>::template funBody<_i>;
  template <fFun_t _f>
  using fFunBody = typename fFunArea_AOS<const Var_t *, Fitness_t *>::template funBody<_f>;

  static inline void defaultInitializeFunctionThatShouldNotBeCalled(Var_t *) {
    constexpr bool
        You_think_you_called_the_default_initialize_function_of_AOS_but_it_is_reserved_as_a_symbol =
            false;
    assert(
        You_think_you_called_the_default_initialize_function_of_AOS_but_it_is_reserved_as_a_symbol);
  }
  static constexpr bool hasParameters = false;
};

#define HEU_MAKE_AOSPARAMETERPACK_TYPES(Base_t) \
  using Args_t = typename Base_t ::Args_t;      \
  using iFun_t = typename Base_t ::iFun_t;      \
  using fFun_t = typename Base_t ::fFun_t;

}  // namespace internal

}  //  namespace heu

#endif  //  HEU_AOSPARAMETERPACK_HPP