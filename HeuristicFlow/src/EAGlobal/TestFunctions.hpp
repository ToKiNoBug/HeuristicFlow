#ifndef HEU_TESTFUNCTIONS_HPP
#define HEU_TESTFUNCTIONS_HPP

#include "InternalHeaderCheck.h"

#include <type_traits>
#include <HeuristicFlow/Global>

namespace heu {
namespace internal {
template <typename Var_t, class Fitness_t = double, class Arg_t = void>
struct SOFunctions {
  static_assert(std::is_floating_point<Fitness_t>::value,
                "Fitness_t must be a floating point number");

  inline static void rastrigin(const Var_t *x, const Arg_t *, Fitness_t *fitness) {
    SOFunctions<Var_t, Fitness_t, void>::rastrigin(x, fitness);
  }

  template <typename unused = void>
  inline static void ackley(const Var_t *x, const Arg_t *, Fitness_t *f) {
    SOFunctions<Var_t>::template ackley<unused>(x, f);
  }
};

template <typename Var_t, class Fitness_t>
struct SOFunctions<Var_t, Fitness_t, void> {
  static_assert(std::is_floating_point<Fitness_t>::value,
                "Fitness_t must be a floating point number");
  inline static void rastrigin(const Var_t *x, Fitness_t *fitness) {
    *fitness = 10 * x->size();
    for (auto xi : *x) {
      *fitness += square(xi) - 10 * std::cos(M_PI * 2 * xi);
    }
  }

  template <typename unused = void>
  inline static void ackley(const Var_t *_x, Fitness_t *f) {
    static_assert((!array_traits<Var_t>::isFixedSize) || (array_traits<Var_t>::sizeCT == 2),
                  "Ackley function requires 2-d array");

    if constexpr (!array_traits<Var_t>::isFixedSize) {
      assert(_x->size() == 2);
    }

    double x = _x->operator[](0), y = _x->operator[](1);
    *f = -20 * std::exp(-0.2 * std::sqrt(0.5 * (x * x + y * y))) -
         std::exp(0.5 * (std::cos(M_2_PI * x) + std::cos(M_2_PI * y))) + 20 + M_E;
  }
};

}  //  namespace internal

template <typename Var_t, class Fitness_t = double, class Arg_t = void>
struct testFunctions
    : public std::conditional<std::is_floating_point<Fitness_t>::value,
                              internal::SOFunctions<Var_t, Fitness_t, Arg_t>, void>::type {};

}  //  namespace heu
#endif
