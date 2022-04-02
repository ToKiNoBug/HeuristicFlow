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

#ifndef HEU_PSOOPTION_HPP
#define HEU_PSOOPTION_HPP

#include <stdint.h>
#include "InternalHeaderCheck.h"
#include <HeuristicFlow/Global>

namespace heu {

// Options to PSO

/**
 * \ingroup HEU_PSO
 * \struct PSOOption
 * \brief Basical option of PSO
 *
 */
struct PSOOption {
 public:
  PSOOption() {
    populationSize = 200;
    maxGeneration = 300;
    maxFailTimes = 100;
    inertiaFactor = 0.8;
    learnFactorP = 2;
    learnFactorG = 2;
  }
  /// size of population, default value is 200
  size_t populationSize;
  /// maximum allowed generation, default value is 300
  size_t maxGeneration;
  /// maximun allowed failtimes, default value is 100
  size_t maxFailTimes;
  /// speed factor, default value is 0.8
  double inertiaFactor;
  /// pBest factor, default value is 2
  double learnFactorP;
  /// gBest factor, default value is 2
  double learnFactorG;
};

}  //  namespace heu

#endif  // HEU_PSOOPTION_HPP
