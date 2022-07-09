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

#ifndef HEU_SOFUNCTIONS_HPP
#define HEU_SOFUNCTIONS_HPP

#include <HeuristicFlow/Global>

#include "testFunctionsCommon.hpp"

#include "InternalHeaderCheck.h"

namespace heu {
namespace internal {

template <typename Var_t, class Fitness_t = double, class Arg_t = void>
struct SOFunctions2 {
  HEU_REPEAT_FUNCTIONS(SOFunctions2, ackley)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, beale)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, GoldSteinPrice)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, booth)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, bukin)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, matyas)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, levy)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, himmelblau)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, easom)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, crossInTray)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, eggHolder)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, holderTable)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, McCormick)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, schaffer2)
  HEU_REPEAT_FUNCTIONS(SOFunctions2, schaffer4)
};

template <typename Var_t, class Fitness_t>
struct SOFunctions2<Var_t, Fitness_t, void> {
 private:
  static_assert(std::is_floating_point<Fitness_t>::value,
                "Fitness_t must be a floating point number");
  static_assert((!array_traits<Var_t>::isFixedSize) || (array_traits<Var_t>::sizeCT == 2),
                "2 testing functions require a 2-d array");
  inline static void assert4Size(const Var_t *_x) noexcept {
    if constexpr (!array_traits<Var_t>::isFixedSize) {
      assert(_x->size() == 2);
    }
  }

 public:
  inline static void ackley(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    auto x = _x->operator[](0), y = _x->operator[](1);
    *f = -20 * std::exp(-0.2 * std::sqrt(0.5 * (x * x + y * y))) -
         std::exp(0.5 * (std::cos(M_2_PI * x) + std::cos(M_2_PI * y))) + 20 + M_E;
  }

  inline static void beale(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = square(1.5 - x + x * y) + square(2.25 - x + x * y * y) + square(2.625 - x + x * y * y * y);
  }

  inline static void GoldsteinPrice(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    const auto part1 =
        1 + square(1 + x + y) * (19 - 14 * x + 3 * x * x - 14 * y + 6 * x * y + 3 * y * y);
    const auto part2 =
        30 + square(2 * x - 3 * y) * (18 - 32 * x + 12 * x + 47 * y - 36 * x * y + 27 * y * y);
    *f = part1 * part2;
  }

  inline static void booth(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = square(x + 2 * y - 7) + square(2 * x + y - 5);
  }

  inline static void bukin(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = 100 * std::sqrt(std::abs(y - 0.01 * x * x)) + 0.01 * std::abs(x + 10);
  }

  inline static void matyas(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = 0.26 * (x * x + y * y) - 0.48 * x * y;
  }

  inline static void levy(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = square(std::sin(3 * M_PI * x)) + square(x - 1) * (1 + square(std::sin(3 * M_PI * y))) +
         square(y - 1) * (1 + square(std::sin(2 * M_PI * y)));
  }

  inline static void himmelblau(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = square(x * x + y - 11) + square(x + y * y - 7);
  }

  inline static void threeHumpCamel(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = 2 * x * x - 1.05 * std::pow(x, 4) + std::pow(x, 6) / 6 + x * y + y * y;
  }

  inline static void easom(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = -std::cos(x) * std::cos(y) * std::exp(-square(x - M_PI) - square(y - M_PI));
  }

  inline static void crossInTray(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    *f = -0.0001 * std::pow(1 + std::abs(std::sin(x) * std::sin(y) *
                                         std::exp(

                                             std::abs(100 - std::sqrt(x * x + y * y) / M_PI)

                                                 )

                                             ),
                            0.1);
  }

  inline static void eggHolder(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);
    const auto y_plus_47 = y + 47;

    *f = -y_plus_47 * std::sin(std::sqrt(std::abs(x / 2 + y_plus_47))) -
         x * std::sin(std::sqrt(std::abs(x - y_plus_47)));
  }

  inline static void holderTable(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);

    *f = -std::abs(std::sin(x) * std::cos(y) *
                   std::exp(

                       std::abs(1 - std::sqrt(x * x - y * y) / M_PI)

                           ));
  }

  inline static void McCormick(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);

    *f = std::sin(x + y) + square(x - y) - 1.5 * x + 2.5 * y + 1;
  }

  inline static void schaffer2(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);

    *f = 0.5 + (square(std::sin(x * x - y * y) - 0.5)) / square(1 + 0.001 * (x * x + y * y));
  }

  inline static void schaffer4(const Var_t *_x, Fitness_t *f) noexcept {
    assert4Size(_x);
    const auto x = _x->operator[](0);
    const auto y = _x->operator[](1);

    *f = 0.5 + (square(std::cos(std::sin(std::abs(x * x - y * y)))) - 0.5) /
                   square(1 + 0.001 * (x * x + y * y));
  }
};

template <typename Var_t, class Fitness_t = double, class Arg_t = void>
struct SOFunctionsX {
  static_assert(std::is_floating_point<Fitness_t>::value,
                "Fitness_t must be a floating point number");

  HEU_REPEAT_FUNCTIONS(SOFunctionsX, rastrigin)
  HEU_REPEAT_FUNCTIONS(SOFunctionsX, sphere)
  HEU_REPEAT_FUNCTIONS(SOFunctionsX, rosenbrock)
  HEU_REPEAT_FUNCTIONS(SOFunctionsX, styblinskiTang)
};

template <typename Var_t, class Fitness_t>
struct SOFunctionsX<Var_t, Fitness_t, void> {
  static_assert(std::is_floating_point<Fitness_t>::value,
                "Fitness_t must be a floating point number");
  inline static void rastrigin(const Var_t *x, Fitness_t *fitness) noexcept {
    *fitness = 10.0 * x->size();
    for (auto xi : *x) {
      *fitness += square(xi) - 10 * std::cos(M_PI * 2 * xi);
    }
  }

  inline static void sphere(const Var_t *x, Fitness_t *f) noexcept {
    if constexpr (array_traits<Var_t>::isEigenClass) {
      *f = (x->array()).square().sum();
    } else {
      double res = 0;
      for (auto i : *x) {
        res += i * i;
      }

      *f = res;
    }
  }

  inline static void rosenbrock(const Var_t *x, Fitness_t *f) noexcept {
    if constexpr (array_traits<Var_t>::isEigenClass) {
      auto square_1_minus_x_top = (1 - x->topRows(x->size() - 1)).square().sum();
      auto minusSquare = (x->topRows(x->size() - 1) - x->bottomRows(x->size() - 1)).square().sum();
      *f = 100 * minusSquare + square_1_minus_x_top;
    } else {
      double val = 0;
      for (int idx = 0; idx < x->size() - 1; idx++) {
        val += 100 * square(x->operator[](idx + 1) - x->operator[](idx)) +
               square(1 - x->operator[](idx));
      }
      *f = val;
    }
  }

  inline static void styblinskiTang(const Var_t *x, Fitness_t *f) noexcept {
    if constexpr (array_traits<Var_t>::isEigenClass) {
      auto arr = x->array();
      *f = (arr.pow(4) - 16 * arr.square() + 5 * arr).sum();
    } else {
      double val = 0;
      for (auto xi : *x) {
        val += xi * xi * xi * xi - 16 * xi * xi + 5 * xi;
      }
      *f = val / 2;
    }
  }
};

}  // namespace internal
}  // namespace heu

#endif  //  HEU_SOFUNCTIONS_HPP