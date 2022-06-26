#ifndef HEU_MOFUNCTIONS_HPP
#define HEU_MOFUNCTIONS_HPP

#include <HeuristicFlow/Global>

#include "testFunctionsCommon.hpp"

#include "InternalHeaderCheck.h"

namespace heu {
namespace internal {

template <typename Var_t, class Fitness_t, class Arg_t = void>
struct MOFunctions12 {
  HEU_REPEAT_FUNCTIONS(MOFunctions12, Schaffer1)
  HEU_REPEAT_FUNCTIONS(MOFunctions12, Schaffer2)
};

template <typename Var_t, class Fitness_t>
struct MOFunctions12<Var_t, Fitness_t, void> {
 private:
  static_assert(sizeMayMatch<Var_t, 1>::value, "MOFunctions12 requires Var_t to be a 1-d array");
  static_assert(sizeMayMatch<Fitness_t, 2>::value,
                "MOFunctions12 requires Fitness_t to be a 2-d array");

  static inline void assert4Size(const Var_t *x) {
    if constexpr (!array_traits<Var_t>::isSizeFixed) {
      assert(x->size() == 1);
    }
  }

  static inline void initializeFitness(Fitness_t *f) {
    if constexpr (!array_traits<Fitness_t>::isSizeFixed) {
      if constexpr (array_traits<Fitness_t>::isEigenClass) {
        f->resize(2, 1);
      } else {
        f->resize(2);
      }
    }
  }

 public:
  static inline void Schaffer1(const Var_t *_x, Fitness_t *f) {
    assert4Size(_x);
    initializeFitness(f);

    const auto x = _x->operator[](0);
    f->operator[](0) = x * x;
    f->operator[](1) = square(x - 2);
  }

  static inline void Schaffer2(const Var_t *_x, Fitness_t *f) {
    assert4Size(_x);
    initializeFitness(f);

    const auto x = _x->operator[](0);
    f->operator[](1) = square(x - 5);

    if (x <= 1) {
      f->operator[](0) = -x;
    } else if (x <= -3) {
      f->operator[](0) = x - 2;
    } else if (x <= 4) {
      f->operator[](0) = 4 - x;
    } else {
      f->operator[](0) = 4 - x;
    }
  }
};

template <typename Var_t, class Fitness_t, class Arg_t = void>
struct MOFunctions22 {
  HEU_REPEAT_FUNCTIONS(MOFunctions22, BinhKorn)
  HEU_REPEAT_FUNCTIONS(MOFunctions22, ChangkongHaimes)
};

template <typename Var_t, class Fitness_t>
struct MOFunctions22<Var_t, Fitness_t, void> {
 private:
  static_assert((!array_traits<Var_t>::isFixedSize) || (array_traits<Var_t>::sizeCT == 2));
  static_assert((!array_traits<Fitness_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT == 2));
  static inline void assert4Size(const Var_t *_x) {
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      assert(_x->size() == 2);
    }
  }

  static inline void initializeFitnessSize(Fitness_t *f) {
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      if constexpr (array_traits<Var_t>::isEigenClass) {
        f->resize(2, 1);
      } else {
        f->resize(2);
      }
    }
  }

 public:
  static inline void BinhKorn(const Var_t *_x, Fitness_t *f) {
    assert4Size(_x);
    initializeFitnessSize(f);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);

    f->operator[](0) = 4 * x * x + 4 * y * y;
    f->operator[](1) = square(x - 5) + square(y - 5);

    const auto penltyG1 = square(x - 5) + y * y - 25;
    const auto penltyG2 = 7.7 - square(x - 8) - square(y + 3);

    if (penltyG1 > 0) f->operator[](0) += 1e5 + penltyG1;

    if (penltyG2 > 0) f->operator[](1) += 1e5 + penltyG2;
  }

  static inline void ChangkongHaimes(const Var_t *_x, Fitness_t *f) {
    assert4Size(_x);
    initializeFitnessSize(f);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);

    f->operator[](0) = 2 + square(x - 2) + square(y - 1);
    f->operator[](1) = 9 * x - square(y - 1);

    const auto penltyG1 = x * x + y * y - 225;

    const auto penltyG2 = x - 3 * y + 10;

    if (penltyG1 > 0) f->operator[](0) += 1e5 + penltyG1;

    if (penltyG2 > 0) f->operator[](1) += 1e5 + penltyG2;
  }

  static inline void Poloni(const Var_t *_x, Fitness_t *f) {
    assert4Size(_x);
    initializeFitnessSize(f);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
#warning Poloni function hasn't been finished yet
    // A1 = 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2)
    constexpr double A1 =
        0.87364856231406409441675015609208455034347519302074844792508553304659184302050078940737876109778881072998046875000000000000000000;
    // A2 = 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2)
    constexpr double A2 =
        2.74857244326863962686864072157634249269337869363896678191405505936963961366448216949720517732203006744384765625000000000000000000
  }
};

template <typename Var_t, class Fitness_t, class Arg_t = void>
struct MOFunctionsX2 {
  HEU_REPEAT_FUNCTIONS(MOFunctionsX2, FonsecaFleming)
  HEU_REPEAT_FUNCTIONS(MOFunctionsX2, Kursawe)
};

template <typename Var_t, class Fitness_t>
struct MOFunctionsX2<Var_t, Fitness_t, void> {
 private:
  static_assert((!array_traits<Fitness_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT == 2));
  static_assert((!array_traits<Var_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT >= 1));

  static inline void assert4Size(const Fitness_t *_x) {
    if constexpr (!array_traits<Fitness_t>::isFixedSize) {
      assert(f->size() == 2);
    }
  }

 public:
  static inline void FonsecaFleming(const Var_t *_x, Fitness_t *f) {
    assert4Size(f);
    const int N = f->size();  // If the size of f is determined at compile-time, the compiler will
                              // be able to optimize it; otherwise N is a read only integer.
    typename array_traits<Fitness_t>::Scalar_t exp1 = 0, exp2 = 0;
    constexpr typename array_traits<Fitness_t>::Scalar_t rSqrtN = 1.0 / std::sqrt(double(N));
    if constexpr (array_traits<Var_t>::isEigenClass) {
      exp1 = -1 * (_x->array() - rSqrtN).square().sum();
      exp2 = -1 * (_x->array() + rSqrtN).square().sum();
    } else {
      for (int idx = 0; idx < N; idx++) {
        exp1 -= sqrt(_x->operator[](idx) - rSqrtN);
        exp2 -= sqrt(_x->operator[](idx) + rSqrtN);
      }
    }

    f->operator[](0) = 1 - std::exp(exp1);
    f->operator[](1) = 1 - std::exp(exp2);
  }

  static inline void Kursawe(const Var_t *_x, Fitness_t *f) {
    assert4Size(f);
    const int N = f->size();
    typename array_traits<Var_t>::Scalar_t f1 = 0, f2 = 0;
    if constexpr (array_traits<Var_t>::isEigenClass) {
      if constexpr (array_traits<Var_t>::isFixedSize) {
        auto topN_sub_1Rows = _x->template topRows<N - 1>();
        auto bottomN_sub_1Rows = _x->template bottomRows<N - 1>();
        f1 = -10 * (-0.2 *

                    (topN_sub_1Rows.array().square() + bottomN_sub_1Rows.array().square())

                        .sqrt())
                       .exp()
                       .sum();

      } else {
        auto topN_sub_1Rows = _x->template topRows(N - 1);
        auto bottomN_sub_1Rows = _x->template bottomRows(N - 1);
        f1 = -10 * (-0.2 *

                    (topN_sub_1Rows.array().square() + bottomN_sub_1Rows.array().square())

                        .sqrt())
                       .exp()
                       .sum();
      }

      f2 = (_x->array().abs().pow(0.8)

            + 5 * (_x->array().cube()).sin()

                )
               .sum();
    } else {
      for (int i = 0; i < N - 1; i++) {
        f1 -= 10 *
              std::exp(-0.2 * std::sqrt(square(_x->operator[](i)) + square(_x->operator[](i + 1))));
        f2 += std::pow(std::abs(_x->operator[](i)), 0.8) +
              5 * std::sin(_x->operator[](i) * _x->operator[](i) * _x->operator[](i));
      }
      const auto lastVarX = _x->operator[](N - 1);
      f2 += std::pow(std::abs(lastVarX), 0.8) + 5 * std::sin(lastVarX * lastVarX * lastVarX);
    }

    f->operator[](0) = f1;
    f->operator[](1) = f2;
  }
};

template <typename Var_t, class Fitness_t, class Arg_t = void>
struct MOFunctionsXX {};

template <typename Var_t, class Fitness_t>
struct MOFunctionsXX<Var_t, Fitness_t, void> {};

}  // namespace internal
}  // namespace heu
#endif  // HEU_MOFUNCTIONS_HPP