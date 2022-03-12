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

#ifndef Heu_TENTCHAOS_H
#define Heu_TENTCHAOS_H

#include <stdint.h>
#include "ChaosBase.hpp"

namespace Heu {

template <uint8_t BitN>
class TentChaosI : public ChaoseBase<uint64_t>
{
public:
    TentChaosI() {
        x=1;
        k=1;
    }

    TentChaosI(uint64_t seed) {
        x=seed&maxVal;
        k=1;
    }

    uint64_t operator()() {
        k<<=2;
        uint64_t g=(x+k)&maxVal;
        if((g>>1)<maxVal) {
            x=(g<<1)+1;
        } else {
            x=(maxVal-g)<<1;
        }
        return x;
    }

    static constexpr uint64_t min() {
        return 0;
    }

    static constexpr uint64_t max() {
        return maxVal;
    }

private:
    uint64_t x;
    uint64_t k;
    constexpr static const uint64_t maxVal=~(~uint64_t(0)<<BitN);

    static_assert (BitN<=64,"BitN mustn't be greater than 64");
    static_assert (BitN>=1,"BitN mustn't be less than 1");
};

}

#endif // TENTCHAOS_H
