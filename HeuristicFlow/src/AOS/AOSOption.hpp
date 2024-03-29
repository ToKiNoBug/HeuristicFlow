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

#ifndef HEU_AOSOPTION_HPP
#define HEU_AOSOPTION_HPP

#include <stdint.h>
#include <stddef.h>

#include "InternalHeaderCheck.h"

namespace heu {

struct AOSOption {
  AOSOption() {
    electronNum = 50;
    maxGeneration = 30;
    maxEarlyStop = 20;
    photonRate = 0.1;
    maxLayerNum = 5;
  }

  ~AOSOption() = default;

  size_t electronNum;
  size_t maxGeneration;
  size_t maxEarlyStop;
  size_t maxLayerNum;
  double photonRate;
};

}  // namespace heu

#endif  //  HEU_AOSOPTION_HPP