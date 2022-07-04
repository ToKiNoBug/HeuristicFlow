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

#ifndef HEU_CONSTANTS_HPP
#define HEU_CONSTANTS_HPP

#include <stdint.h>
#include <limits>

#include "InternalHeaderCheck.h"

namespace heu {

namespace internal {

/**
 * \ingroup HEU_GLOBAL
 * \brief Positive infinity float
 *
 */
constexpr float pinfF = std::numeric_limits<float>::infinity();

/**
 * \ingroup HEU_GLOBAL
 * \brief Positive infinity double
 *
 */
constexpr double pinfD = std::numeric_limits<double>::infinity();

/**
 * \ingroup HEU_GLOBAL
 * \brief Negativev infinity float
 *
 */
constexpr float ninfF = -pinfF;

/**
 * \ingroup HEU_GLOBAL
 * \brief Negative infinity double
 *
 */
constexpr double ninfD = -pinfD;

}  //  namespace internal

}  //  namespace heu

#endif  // HEU_CONSTANTS_HPP
