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

#ifndef HEU_FUNBODY_HPP
#define HEU_FUNBODY_HPP

#include <type_traits>

#include "InternalHeaderCheck.h"

namespace heu {

/**
 * \ingroup HEU_GLOBAL
 * \brief This marco metafunction is used to generate a struct that adpatively mantain a function
 * pointer.
 *
 *
 * If a function name is provided at compile time, it can execute the function,
 * while if nullptr is assigned, it can store a function pointer and get its value
 * at runtime.
 *
 */
#define HEU_MAKE_FUNAREA(funFlag, Suffix)                                                         \
  template <typename... a>                                                                        \
  class funFlag##Area_##Suffix {                                                                  \
   public:                                                                                        \
    using funPtr_t = void (*)(a...);                                                              \
                                                                                                  \
   private:                                                                                       \
    template <void (*_fun)(a...)>                                                                 \
    class funBodyCT {                                                                             \
     public:                                                                                      \
      inline constexpr funPtr_t funFlag() const { return _fun; }                                  \
      inline void run##funFlag(a... _a) const { _fun(_a...); }                                    \
      constexpr funPtr_t funFlag##AtCompileTime = _fun;                                           \
                                                                                                  \
     private:                                                                                     \
      static_assert(_fun != nullptr, "Template function mustn't be nullptr");                     \
    };                                                                                            \
                                                                                                  \
    class funBodyRT {                                                                             \
     public:                                                                                      \
      funBodyRT() { _funPtr = nullptr; }                                                          \
      inline funPtr_t funFlag() const { return _funPtr; }                                         \
      inline void run##funFlag(a... _a) const { _funPtr(_a...); }                                 \
      inline void set##funFlag(funPtr_t __) { _funPtr = __; }                                     \
      constexpr funPtr_t funFlag##AtCompileTime = nullptr;                                        \
                                                                                                  \
     protected:                                                                                   \
      funPtr_t _funPtr;                                                                           \
    };                                                                                            \
                                                                                                  \
   public:                                                                                        \
    template <void (*_fun)(a...)>                                                                 \
    using funBody = typename std::conditional<_fun == nullptr, funBodyRT, funBodyCT<_fun>>::type; \
  };

}  //  namespace heu

#endif  // HEU_FUNBODY_HPP
