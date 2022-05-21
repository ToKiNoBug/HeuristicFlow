

#ifndef HEU_FUNBODY_HPP
#define HEU_FUNBODY_HPP

#include <type_traits>

#include "InternalHeaderCheck.h"

namespace heu {

/**
 * \ingroup HEU_GLOBAL
 * \brief This marco metafunction is used to generate a struct that adpatively mantain a function pointer.
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
