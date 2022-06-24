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

#ifndef HEU_AOSDEFAULTELECTRON_HPP
#define HEU_AOSDEFAULTELECTRON_HPP

#include "InternalHeaderCheck.hpp"

namespace heu {
namespace internal {
template <typename Var_t, class Fitness_t>
struct DefaultElectron {
  DefaultElectron() { isComputed = false; }
  ~DefaultElectron() = default;
  Var_t state;
  Fitness_t energy;
  bool isComputed;
  inline void setUncomputed() noexcept { isComputed = false; }
};
}  // namespace internal
}  // namespace heu

#endif  //  HEU_AOSDEFAULTELECTRON_HPP