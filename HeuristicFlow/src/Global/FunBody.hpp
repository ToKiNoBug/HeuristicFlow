// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_FUNBODY_HPP
#define EIGEN_HEU_FUNBODY_HPP

#include <type_traits>

#include "InternalHeaderCheck.h"

namespace Eigen {

/**
 * \brief This marco metafunction is used to generate a struct
 * that adpatively mantain a function pointer. If a function name
 * is provided at compile time, it can execute it, while if nullptr
 * is assigned, it can store a function pointer and get its value
 * at runtime.
 *
 */
#define EIGEN_HEU_MAKE_FUNAREA(funFlag, FunFlag, Suffix)                                          \
  template <typename return_t = void, typename... a>                                              \
  class funFlag##Area_##Suffix {                                                                  \
   public:                                                                                        \
    using funPtr_t = return_t (*)(a...);                                                          \
                                                                                                  \
   private:                                                                                       \
    template <return_t (*_fun)(a...)>                                                             \
    class funBodyCT {                                                                             \
     public:                                                                                      \
      inline constexpr funPtr_t funFlag() const { return _fun; }                                  \
      inline return_t run##FunFlag(a... _a) const { return _fun(_a...); }                         \
                                                                                                  \
     private:                                                                                     \
      static_assert(_fun != nullptr, "Template function mustn't be nullptr");                     \
    };                                                                                            \
                                                                                                  \
    class funBodyRT {                                                                             \
     public:                                                                                      \
      funBodyRT() { _funPtr = nullptr; }                                                          \
      inline funPtr_t funFlag() const { return _funPtr; }                                         \
      inline return_t run##FunFlag(a... _a) const { return _funPtr(_a...); }                      \
      inline void set##FunFlag(funPtr_t __) { _funPtr = __; }                                     \
                                                                                                  \
     protected:                                                                                   \
      funPtr_t _funPtr;                                                                           \
    };                                                                                            \
                                                                                                  \
   public:                                                                                        \
    template <return_t (*_fun)(a...)>                                                             \
    using funBody = typename std::conditional<_fun == nullptr, funBodyRT, funBodyCT<_fun>>::type; \
  };

}  //  namespace Eigen

#endif  // EIGEN_HEU_FUNBODY_HPP
