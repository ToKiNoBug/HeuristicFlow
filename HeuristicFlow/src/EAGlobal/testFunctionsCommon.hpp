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

#ifndef HEU_TESTFUNCTIONSCOMMON_HPP
#define HEU_TESTFUNCTIONSCOMMON_HPP
namespace heu {
namespace internal {

#define HEU_REPEAT_FUNCTIONS(className, functionName)                                     \
  inline static void functionName(const Var_t *x, const Arg_t *, Fitness_t *f) noexcept { \
    className<Var_t, Fitness_t, void>::functionName(x, f);                                \
  }

}  // namespace internal
}  // namespace heu
#endif  // HEU_TESTFUNCTIONSCOMMON_HPP