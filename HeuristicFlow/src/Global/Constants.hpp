/*
 Copyright © 2022  TokiNoBug
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

#ifndef Heu_CONSTANTS_HPP
#define Heu_CONSTANTS_HPP

#include <stdint.h>
#include <limits>

namespace Heu {

///Size identifier for dynamic size (fitness or var)
const size_t Runtime = 0;

///inf value for float
const float pinfF=std::numeric_limits<float>::infinity();

///inf value for double
const double pinfD=std::numeric_limits<double>::infinity();

///negative inf value for float
const float nInfF=-pinfF;

///negative inf value for double
const double ninfD=-pinfD;

}

#endif // Heu_CONSTANTS_HPP
