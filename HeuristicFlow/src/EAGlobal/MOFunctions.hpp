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

#ifndef HEU_MOFUNCTIONS_HPP
#define HEU_MOFUNCTIONS_HPP

#include <HeuristicFlow/Global>

#include "testFunctionsCommon.hpp"

#include "InternalHeaderCheck.h"

#include <iostream>

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

  static inline void assert4Size(const Var_t *x) noexcept {
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      assert(x->size() == 1);
    }
  }

  static inline void initializeFitness(Fitness_t *f) noexcept {
    if constexpr (!array_traits<Fitness_t>::isFixedSize) {
      if constexpr (array_traits<Fitness_t>::isEigenClass) {
        f->resize(2, 1);
      } else {
        f->resize(2);
      }
    }
  }

 public:
  static inline void Schaffer1(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    initializeFitness(f);

    const auto x = _x->operator[](0);
    f->operator[](0) = x * x;
    f->operator[](1) = square(x - 2);
  }

  static inline void Schaffer2(const Var_t *_x, Fitness_t *f) noexcept {
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
  HEU_REPEAT_FUNCTIONS(MOFunctions22, Poloni)
};

template <typename Var_t, class Fitness_t>
struct MOFunctions22<Var_t, Fitness_t, void> {
 private:
  static_assert((!array_traits<Var_t>::isFixedSize) || (array_traits<Var_t>::sizeCT == 2));
  static_assert((!array_traits<Fitness_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT == 2));
  static inline void assert4Size(const Var_t *_x) noexcept {
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      assert(_x->size() == 2);
    }
  }

  static inline void initializeFitnessSize(Fitness_t *f) noexcept {
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      if constexpr (array_traits<Var_t>::isEigenClass) {
        f->resize(2, 1);
      } else {
        f->resize(2);
      }
    }
  }

 public:
  static inline void BinhKorn(const Var_t *_x, Fitness_t *f) noexcept {
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

  static inline void ChangkongHaimes(const Var_t *_x, Fitness_t *f) noexcept {
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

  static inline void Poloni(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    initializeFitnessSize(f);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    // A1 = 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2)
    constexpr double A1 =
        0.87364856231406409441675015609208455034347519302074844792508553304659184302050078940737876109778881072998046875000000000000000000;
    // A2 = 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2)
    constexpr double A2 =
        2.74857244326863962686864072157634249269337869363896678191405505936963961366448216949720517732203006744384765625000000000000000000;

    const auto B1 = 0.5 * std::sin(x) - 2 * std::cos(x) + std::sin(y) - 1.5 * std::cos(y);
    const auto B2 = 1.5 * std::sin(x) - std::cos(x) + 2 * std::sin(y) - 0.5 * std::cos(y);

    f->operator[](0) = 1 + square(A1 - B1) + square(A2 - B2);
    f->operator[](1) = square(x + 3) + square(y + 1);
  }
};

template <typename Var_t, class Fitness_t, class Arg_t = void>
struct MOFunctions23 {
  HEU_REPEAT_FUNCTIONS(MOFunctions23, viennet)
};

template <typename Var_t, class Fitness_t>
struct MOFunctions23<Var_t, Fitness_t, void> {
 private:
  static_assert((!array_traits<Var_t>::isFixedSize) || (array_traits<Var_t>::sizeCT == 2));
  static_assert((!array_traits<Fitness_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT == 3));
  static inline void assert4Size(const Var_t *_x) noexcept {
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      assert(_x->size() == 2);
    }
  }

  static inline void initializeFitnessSize(Fitness_t *f) noexcept {
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      if constexpr (array_traits<Var_t>::isEigenClass) {
        f->resize(3, 1);
      } else {
        f->resize(3);
      }
    }
  }

 public:
  static inline void viennet(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    initializeFitnessSize(f);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);

    const auto x2_plus_y2 = x * x + y * y;
    f->operator[](0) = x2_plus_y2 / 2 + std::sin(x2_plus_y2);
    f->operator[](1) = square(2 * x - 2 * y + 4) / 8 + square(x - y + 1) / 27 + 15;
    f->operator[](2) = 1 / (x2_plus_y2 + 1) - 1.1 * std::exp(-x2_plus_y2);
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

  static inline void assert4Size(const Fitness_t *f) noexcept {
    if constexpr (!array_traits<Fitness_t>::isFixedSize) {
      assert(f->size() == 2);
    }
  }

 public:
  static inline void FonsecaFleming(const Var_t *_x, Fitness_t *f) noexcept {
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

  static inline void Kursawe(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(f);
    const int N = int(f->size());
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
struct DTLZ1to7 {
  HEU_REPEAT_FUNCTIONS(DTLZ1to7, DTLZ1)
  HEU_REPEAT_FUNCTIONS(DTLZ1to7, DTLZ2)
  HEU_REPEAT_FUNCTIONS(DTLZ1to7, DTLZ3)

  template <int alpha = 100>
  inline static void DTLZ4(const Var_t *x, const Arg_t *, Fitness_t *f) noexcept {
    DTLZ1to7<Var_t, Fitness_t, void>::template DTLZ4<alpha>(x, f);
  }

  HEU_REPEAT_FUNCTIONS(DTLZ1to7, DTLZ5)
  HEU_REPEAT_FUNCTIONS(DTLZ1to7, DTLZ6)
  HEU_REPEAT_FUNCTIONS(DTLZ1to7, DTLZ7)
};

template <typename Var_t, class Fitness_t>
struct sizeMayMatchDTLZ1to7 {
  static constexpr bool isAllFixed =
      array_traits<Var_t>::isFixedSize && array_traits<Fitness_t>::isFixedSize;
  static_assert((!array_traits<Fitness_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT));

  static constexpr bool value =
      (isAllFixed) ? (array_traits<Var_t>::sizeCT - array_traits<Fitness_t>::sizeCT + 1 > 0)
                   : (true);
};

template <typename Var_t, class Fitness_t>
struct sizeMayMatchDTLZ89 {
  static constexpr bool isAllFixed =
      array_traits<Var_t>::isFixedSize && array_traits<Fitness_t>::isFixedSize;
  static_assert((!array_traits<Fitness_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT));

  static constexpr bool value =
      (isAllFixed) ? (array_traits<Var_t>::sizeCT > array_traits<Fitness_t>::sizeCT) : (true);
};

template <typename Var_t, class Fitness_t>
struct DTLZ1to7<Var_t, Fitness_t, void> {
 private:
  static_assert((!array_traits<Var_t>::isFixedSize) || (array_traits<Var_t>::sizeCT > 1));
  static_assert((!array_traits<Fitness_t>::isFixedSize) || (array_traits<Fitness_t>::sizeCT > 1));

  inline static void assert4Size(const Var_t *_x, const Fitness_t *f) noexcept {
    // Var_t.size = N,Fitness_t.size =M. k=N-M+1>0
    if constexpr ((array_traits<Var_t>::isFixedSize) && (array_traits<Fitness_t>::isFixedSize)) {
      constexpr int N = array_traits<Var_t>::sizeCT;
      constexpr int M = array_traits<Fitness_t>::sizeCT;
      static_assert(N - M + 1 > 0, "DTLZ1 to DTLZ7 problems require N-M+1 to be positive");
    } else {
      const int N = int(_x->size());
      const int M = int(f->size());
      const bool DTLZ1_to_DTLZ7_problems_require__N_minus_M_plus_1__to_be_positive = N - M + 1 > 0;
      assert(DTLZ1_to_DTLZ7_problems_require__N_minus_M_plus_1__to_be_positive);
    }
  }

  using scalar_t = typename array_traits<Fitness_t>::Scalar_t;
  static constexpr bool isAllFixed =
      array_traits<Var_t>::isFixedSize && array_traits<Fitness_t>::isFixedSize;

  inline static auto computeXM(const Var_t *x, const Fitness_t *f) noexcept {
    if constexpr (isAllFixed) {
      constexpr int N = array_traits<Var_t>::sizeCT;
      constexpr int M = array_traits<Fitness_t>::sizeCT;
      return x->template tail<N - M + 1>();
    } else {
      const int N = int(x->size());
      const int M = int(f->size());
      return x->tail(N - M + 1);
    }
  }

  // Type of vector X except XM. It's of the same vector type with Var_t
  using XPrevPartVec_t = typename std::conditional<
      array_traits<Var_t>::isEigenClass,

      typename std::conditional<
          array_traits<Fitness_t>::isFixedSize,
          Eigen::Array<typename array_traits<Var_t>::Scalar_t, array_traits<Fitness_t>::sizeCT - 1,
                       1>,
          Eigen::Array<typename array_traits<Var_t>::Scalar_t, Eigen::Dynamic, 1>>::type,

      typename std::conditional<
          array_traits<Fitness_t>::isFixedSize,
          std::array<typename array_traits<Var_t>::Scalar_t, array_traits<Fitness_t>::sizeCT - 1>,
          std::vector<typename array_traits<Var_t>::Scalar_t>>::type

      >::type;

 public:
  // x=[x1;x2;...;xM-1;vec(xM)]
  // xM is the last N-M+1 elements in x
  // index of xM starts from M-1 and end with N-1
  inline static void DTLZ1(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);

    scalar_t g = 0;
    const int N = int(x->size());
    const int M = int(f->size());
    if constexpr (array_traits<Var_t>::isEigenClass) {
      auto xm = computeXM(x, f);
      g = 100 * (

                    xm.size() +

                    ((xm - 0.5).square() - (20 * M_PI * (xm - 0.5)).cos()).sum()

                );
    } else {
      g = 100 * (N - M + 1);
      for (int idx = M - 1; idx < N; idx++) {
        const scalar_t xi_sub_half = x->operator[](idx) - 0.5;
        g += 100 * (xi_sub_half * xi_sub_half

                    - std::cos(20 * M_PI * xi_sub_half));
      }
    }

    scalar_t prodCum = 0.5 * (1 + g);

    for (int objIdx = M - 1; objIdx > 0; objIdx--) {
      f->operator[](objIdx) = prodCum * (1 - x->operator[](M - 1 - objIdx));
      prodCum *= x->operator[](M - 1 - objIdx);
    }
    f->operator[](0) = prodCum;
  }

  inline static void DTLZ2(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);

    const int N = int(x->size());
    const int M = int(f->size());

    scalar_t g = 0;

    if constexpr (array_traits<Var_t>::isEigenClass) {
      auto xm = computeXM(x, f);
      g = (xm - 0.5).square().sum();
    } else {
      for (int idx = M - 1; idx < N; idx++) {
        g += square(x->operator[](idx) - 0.5);
      }
    }

    scalar_t prodCum = 1 + g;
    for (int objIdx = M - 1; objIdx > 0; objIdx--) {
      f->operator[](objIdx) = prodCum * std::sin(M_PI / 2 * x->operator[](M - 1 - objIdx));
      prodCum *= std::cos(M_PI / 2 * x->operator[](M - 1 - objIdx));
    }
    f->operator[](0) = prodCum;
  }

  inline static void DTLZ3(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);

    const int N = int(x->size());
    const int M = int(f->size());

    scalar_t g = 0;

    if constexpr (array_traits<Var_t>::isEigenClass) {
      auto xm = computeXM(x, f);
      g = 100 * (

                    xm.size() +

                    ((xm - 0.5).square() - (20 * M_PI * (xm - 0.5)).cos()).sum()

                );
    } else {
      g = 100 * (N - M + 1);
      for (int idx = M - 1; idx < N; idx++) {
        const scalar_t xi_sub_half = x->operator[](idx) - 0.5;
        g += 100 * (xi_sub_half * xi_sub_half

                    - std::cos(20 * M_PI * xi_sub_half));
      }
    }

    scalar_t prodCum = 1 + g;
    for (int objIdx = M - 1; objIdx > 0; objIdx--) {
      f->operator[](objIdx) = prodCum * std::sin(M_PI / 2 * x->operator[](M - 1 - objIdx));
      prodCum *= std::cos(M_PI / 2 * x->operator[](M - 1 - objIdx));
    }
    f->operator[](0) = prodCum;
  }

  template <int alpha = 100>
  inline static void DTLZ4(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);

    const int N = int(x->size());
    const int M = int(f->size());

    scalar_t g = 0;

    if constexpr (array_traits<Var_t>::isEigenClass) {
      auto xm = computeXM(x, f);
      g = (xm - 0.5).square().sum();
    } else {
      for (int idx = M - 1; idx < N; idx++) {
        g += square(x->operator[](idx) - 0.5);
      }
    }

    scalar_t prodCum = 1 + g;
    for (int objIdx = M - 1; objIdx > 0; objIdx--) {
      f->operator[](objIdx) =
          prodCum * std::sin(M_PI / 2 * std::pow(x->operator[](M - 1 - objIdx), alpha));
      prodCum *= std::cos(M_PI / 2 * std::pow(x->operator[](M - 1 - objIdx), alpha));
    }
    f->operator[](0) = prodCum;
  }

  inline static void DTLZ5(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);

    const int N = int(x->size());
    const int M = int(f->size());

    scalar_t gXm = 0, gr = 0;

    if (array_traits<Var_t>::isEigenClass) {
      auto xm = computeXM(x, f);
      gXm = (xm - 0.5).square().sum();
      gr = xm.square().sum();
    } else {
      for (int idx = M - 1; idx < N; idx++) {
        gr += square(x->operator[](idx));
        gXm += square(x->operator[](idx) - 0.5);
      }
    }

    XPrevPartVec_t theta_times_PI_div_2;
    // compute theta;
    if constexpr (array_traits<XPrevPartVec_t>::isEigenClass) {  // Eigen's vector
      if constexpr (array_traits<Fitness_t>::isFixedSize) {
        // theta = (1 + 2 * gr * x->template topRows<M - 1>()) * M_PI / (4 * (1 + gr));
        constexpr int M = array_traits<Fitness_t>::sizeCT;  // reload the M from const to constexpr
        theta_times_PI_div_2 =
            M_PI / 2 * (1 + 2 * gr * x->template head<M - 1>()) * M_PI / (4 * (1 + gr));
      } else {
        // theta = (1 + 2 * gr * x->topRows(M - 1)) * M_PI / (4 * (1 + gr));
        theta_times_PI_div_2 = M_PI / 2 * (1 + 2 * gr * x->head(M - 1)) * M_PI / (4 * (1 + gr));
      }

      // theta[0] = M_PI / 2 * x->operator[](0);
      theta_times_PI_div_2[0] = M_PI / 2 * M_PI / 2 * x->operator[](0);
    } else {
      if constexpr (!array_traits<XPrevPartVec_t>::isFixedSize) {  // std::array or std::vector
        // theta.resize(M - 1);
        theta_times_PI_div_2.resize(M - 1);
        // resize for std::vector
      }

      // theta[0] = M_PI / 2 * x->operator[](0);
      theta_times_PI_div_2[0] = M_PI / 2 * M_PI / 2 * x->operator[](0);

      for (int idx = 1; idx < M - 1; idx++) {
        // theta[idx] = (1 + 2 * gr * x->operator[](idx)) * M_PI / (4 * (1 + gr));
        theta_times_PI_div_2[idx] = M_PI / 2 *

                                    (1 + 2 * gr * x->operator[](idx)) * M_PI / (4 * (1 + gr));
      }
    }

    scalar_t prodCum = 1 + gXm;

    for (int objIdx = M - 1; objIdx > 0; objIdx--) {
      f->operator[](objIdx) = prodCum * std::sin(theta_times_PI_div_2[M - 1 - objIdx]);
      prodCum *= std::cos(theta_times_PI_div_2[M - 1 - objIdx]);
    }
    f->operator[](0) = prodCum;
  }

  inline static void DTLZ6(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);

    const int N = int(x->size());
    const int M = int(f->size());

    scalar_t gXm = 0, gr = 0;

    if (array_traits<Var_t>::isEigenClass) {
      auto xm = computeXM(x, f);
      gr = xm.square().sum();
      gXm = xm.pow(0.1).sum();
    } else {
      for (int idx = M - 1; idx < N; idx++) {
        gr += square(x->operator[](idx));
        gXm += std::pow(x->operator[](idx), 0.1);
      }
    }

    XPrevPartVec_t theta_times_PI_div_2;
    // compute theta;
    if constexpr (array_traits<XPrevPartVec_t>::isEigenClass) {  // Eigen's vector
      if constexpr (array_traits<Fitness_t>::isFixedSize) {
        // theta = (1 + 2 * gr * x->template topRows<M - 1>()) * M_PI / (4 * (1 + gr));
        constexpr int M = array_traits<Fitness_t>::sizeCT;  // reload the M from const to constexpr
        theta_times_PI_div_2 =
            M_PI / 2 * (1 + 2 * gr * x->template head<M - 1>()) * M_PI / (4 * (1 + gr));
      } else {
        // theta = (1 + 2 * gr * x->topRows(M - 1)) * M_PI / (4 * (1 + gr));
        theta_times_PI_div_2 = M_PI / 2 * (1 + 2 * gr * x->head(M - 1)) * M_PI / (4 * (1 + gr));
      }

      // theta[0] = M_PI / 2 * x->operator[](0);
      theta_times_PI_div_2[0] = M_PI / 2 * M_PI / 2 * x->operator[](0);
    } else {
      if constexpr (!array_traits<XPrevPartVec_t>::isFixedSize) {  // std::array or std::vector
        // theta.resize(M - 1);
        theta_times_PI_div_2.resize(M - 1);
        // resize for std::vector
      }

      // theta[0] = M_PI / 2 * x->operator[](0);
      theta_times_PI_div_2[0] = M_PI / 2 * M_PI / 2 * x->operator[](0);

      for (int idx = 1; idx < M - 1; idx++) {
        // theta[idx] = (1 + 2 * gr * x->operator[](idx)) * M_PI / (4 * (1 + gr));
        theta_times_PI_div_2[idx] = M_PI / 2 *

                                    (1 + 2 * gr * x->operator[](idx)) * M_PI / (4 * (1 + gr));
      }
    }

    scalar_t prodCum = 1 + gXm;

    for (int objIdx = M - 1; objIdx > 0; objIdx--) {
      f->operator[](objIdx) = prodCum * std::sin(theta_times_PI_div_2[M - 1 - objIdx]);
      prodCum *= std::cos(theta_times_PI_div_2[M - 1 - objIdx]);
    }
    f->operator[](0) = prodCum;
  }

  inline static void DTLZ7(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);

    const int N = int(x->size());
    const int M = int(f->size());

    // compute the first M-1 objective values
    if constexpr (array_traits<Var_t>::isEigenClass && array_traits<Fitness_t>::isEigenClass) {
      if constexpr (array_traits<Fitness_t>::isFixedSize) {
        // constexpr int N = array_traits<Var_t>::sizeCT;
        constexpr int M = array_traits<Fitness_t>::sizeCT;
        f->template head<M - 1>() = x->template head<M - 1>();
      } else {
        f->head(M - 1) = x->head(M - 1);
      }
    } else {
      for (int idx = 0; idx < M - 2; idx++) {
        f->operator[](idx) = x->operator[](idx);
      }
    }

    scalar_t gXm = 0;
    if constexpr (array_traits<Var_t>::isEigenClass) {
      auto xm = computeXM(x, f);
      gXm = 1 + 9.0 / xm.size() * xm.sum();
    } else {
      // gXm = 1;
      for (int idx = M - 1; idx < N; idx++) {
        gXm += x->operator[](idx);
      }
      gXm *= 9.0 / (N - M + 1);
      gXm += 1;
    }

    scalar_t h = 0;
    for (int objIdx = 0; objIdx < M - 1; objIdx++)
      h -= f->operator[](objIdx) / (1 + gXm) * (1 + std::sin(3 * M_PI * f->operator[](objIdx)));

    h += M;

    f->operator[](M - 1) = (1 + gXm) * h;
  }
};

template <typename Var_t, class Fitness_t, class Arg_t = void>
struct DTLZ89 {
  HEU_REPEAT_FUNCTIONS(DTLZ89, DTLZ8)
  HEU_REPEAT_FUNCTIONS(DTLZ89, DTLZ9)
};

template <typename Var_t, class Fitness_t>
struct DTLZ89<Var_t, Fitness_t, void> {
 private:
  inline static void assert4Size(const Var_t *_x, const Fitness_t *f) noexcept {
    // Var_t.size = N,Fitness_t.size =M. N>M must be true
    if constexpr (isAllFixed) {
      constexpr int N = array_traits<Var_t>::sizeCT;
      constexpr int M = array_traits<Fitness_t>::sizeCT;
      static_assert(N > M,
                    "DTLZ8 and DTLZ9 require the dimensions of determine variable to be greater "
                    "than the number of objectives.");
    } else {
      const int N = int(_x->size());
      const int M = int(f->size());
      const bool
          DTLZ8_and_DTLZ9_require_the_dimensions_of_determine_variable_to_be_greater_than_the_number_of_objectives =
              N > M;

      assert(
          DTLZ8_and_DTLZ9_require_the_dimensions_of_determine_variable_to_be_greater_than_the_number_of_objectives);
    }
  }

  using scalar_t = typename array_traits<Fitness_t>::Scalar_t;
  static constexpr bool isAllFixed =
      array_traits<Var_t>::isFixedSize && array_traits<Fitness_t>::isFixedSize;

  template <int varIdx, int endIdx>
  inline static scalar_t accumulateXi(const Var_t *x) noexcept {
    if constexpr (varIdx < endIdx) {
      return x->operator[](varIdx) + accumulateXi<varIdx + 1, endIdx>(x);
    } else {
      return x->operator[](endIdx);
    }
  }

  inline static scalar_t accumulateXi(const Var_t *x, const int begIdx, const int endIdx) noexcept {
    if constexpr (array_traits<Var_t>::isEigenClass) {
      return x->segment(begIdx, endIdx - begIdx + 1).sum();
    } else {
      scalar_t sum = 0;
      for (int idx = begIdx; idx <= endIdx; idx++) {
        sum += x->operator[](idx);
      }
      return sum;
    }
  }

  template <int objIdx>
  inline static void computeFitness8(const Var_t *x, Fitness_t *f) noexcept {
    static_assert(array_traits<Fitness_t>::isFixedSize);
    constexpr int M = array_traits<Fitness_t>::sizeCT;
    if constexpr (objIdx < M) {
      if constexpr (array_traits<Var_t>::isFixedSize) {
        constexpr int N = array_traits<Var_t>::sizeCT;
        constexpr int begIdx = objIdx * N / M;
        constexpr int tempEndIdx = (objIdx + 1) * N / M;
        constexpr int endIdx = (tempEndIdx >= N) ? (N - 1) : (tempEndIdx);

        f->operator[](objIdx) = accumulateXi<begIdx, endIdx>(x);

      } else {
        const int N = x->size();
        const int begIdx = objIdx * N / M;
        const int endIdx = std::min<int>((objIdx + 1) * N / M, N - 1);
        f->operator[](objIdx) = accumulateXi(x, begIdx, endIdx);
      }
      computeFitness8<objIdx + 1>(x, f);
    } else {
      return;
    }
  }

  inline static void computeFitness8(const Var_t *x, Fitness_t *f) noexcept {
    // Here are two implementations for fixed and dynamic Fitness_t.
    // For compile-time M, expand the 2-rank loop into 1-rank loop or 0-rank loop by templates.
    // Use EIGEN_FORCED_INLINE as many possible.
    // For runtime M, use the 2-rank loop.
    if constexpr (array_traits<Fitness_t>::isFixedSize) {
      computeFitness8<0>(x, f);
    } else {
      const int N = x->size();
      const int M = f->size();
      for (int objIdx = 0; objIdx < M; objIdx++) {
        f->operator[](objIdx) =
            accumulateXi(x, objIdx * N / M, std::min<int>((objIdx + 1) * N / M, N - 1));
      }
    }
  }

  template <int i, int j>  // 0<=i<j<=M-2
  inline static scalar_t findMin2SumForDTLZ8(Fitness_t *f) noexcept {
    static_assert(array_traits<Fitness_t>::isFixedSize);
    constexpr int M = array_traits<Fitness_t>::sizeCT;
    static_assert(i >= 0);
    static_assert(j > i);

    // the outer expansion adjust i from 0 to M-3
    // and the inner expansion adjust j from i+1 to M-2
    if constexpr (i < M - 3) {
      if constexpr (j < M - 1) {
        return std::min(f->operator[](i) + f->operator[](j), findMin2SumForDTLZ8<i, j + 1>(f));
      } else {
        return findMin2SumForDTLZ8<i + 1, i + 2>(f);
      }
    } else {
      static_assert(i == M - 3, "Wrong in programming");
      return f->operator[](i) + f->operator[](j);
      // here i should be M-3
    }
  }

  inline static scalar_t findMin2SumForDTLZ8(Fitness_t *f) noexcept {
    if constexpr (array_traits<Fitness_t>::isFixedSize) {
      return findMin2SumForDTLZ8<0, 1>(f);
    } else {
      const int M = f->size();
      scalar_t result = pinfD;
      for (int i = 0; i < M - 2; i++) {
        for (int j = i + 1; j < M - 1; j++) {
          result = std::min(result, f->operator[](i) + f->operator[](j));
        }
      }

      return result;
    }
  }

 public:
  inline static void DTLZ8(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);
    const int N = int(x->size());
    const int M = int(f->size());
    const int floor_N_div_M = N / M;

    computeFitness8(x, f);

    for (auto &val : *f) {
      val /= floor_N_div_M;
    }

    Fitness_t g;
    if constexpr (!array_traits<Fitness_t>::isFixedSize) {
      if constexpr (array_traits<Fitness_t>::isEigenClass) {
        g.resize(f->size(), 1);  // Eigen::ArrayX
      } else {
        g.resize(f->size());  // std::vector
      }
    }

    for (int objIdx = 0; objIdx < M - 1; objIdx++) {
      g[objIdx] = f->operator[](M - 1) + 4 * f->operator[](objIdx) - 1;
    }

    g[M - 1] = 2 * f->operator[](M - 1) + findMin2SumForDTLZ8(f) - 1;

    for (int objIdx = 0; objIdx < M; objIdx++) {
      if (g[objIdx] < 0) {
        f->operator[](objIdx) += 1e3 + (-g[objIdx]);
      }
    }
  }

  inline static void DTLZ9(const Var_t *x, Fitness_t *f) noexcept {
    assert4Size(x, f);

    const int N = int(x->size());
    const int M = int(f->size());

    for (int objIdx = 0; objIdx < M; objIdx++) {
      f->operator[](objIdx) = 0;
      for (int varIdx = objIdx * N / M; varIdx <= std::min<int>(N, (objIdx + 1) * N / M);
           varIdx++) {
        f->operator[](objIdx) += std::pow(x->operator[](varIdx), 0.1);
      }
    }
    // apply penlties
    for (int objIdx = 0; objIdx < M - 1; objIdx++) {
      auto gOfObjIdx = square(f->operator[](M)) + square(f->operator[](objIdx)) - 1;

      if (gOfObjIdx < 0) {
        f->operator[](objIdx) += 1e3 + (-gOfObjIdx);
      }
    }
  }
};

}  // namespace internal
}  // namespace heu
#endif  // HEU_MOFUNCTIONS_HPP