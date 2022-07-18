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

#ifndef HEU_DEFAULTGENETYPE_HPP
#define HEU_DEFAULTGENETYPE_HPP

#include "InternalHeaderCheck.h"

#include "../../Global"
#include "../../EAGlobal"
#include "GAAbstract.hpp"

namespace heu {

namespace internal {
/**
 * \ingroup HEU_GENETIC
 * \brief Default type of gene for GA
 *
 * \tparam Var_t
 * \tparam Fitness_t
 */
template <typename Var_t, typename Fitness_t>
class DefaultGene_t {
 public:
  using fastFitness_t = typename std::conditional<(sizeof(Fitness_t) > sizeof(double)),
                                                  const Fitness_t &, Fitness_t>::type;

  Var_t self;          ///< Value of decision variable
  Fitness_t _Fitness;  ///< Value of fitness
  bool _isCalculated;  ///< Whether the fitness is computed

  inline bool isCalculated() const noexcept {
    return _isCalculated;
  }  ///< If the fitness is computed
  inline void setUncalculated() noexcept {
    _isCalculated = false;
  }  ///< Set the fitness to be uncomputed
  inline fastFitness_t fitness() const noexcept { return _Fitness; }  ///< Get fitness
};

}  // namespace internal

}  // namespace heu

#endif  //  HEU_DEFAULTGENETYPE_HPP