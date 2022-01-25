/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef TENTCHAOS_H
#define TENTCHAOS_H

#include <stdint.h>

namespace OptimT {

template <uint8_t BitN>
class TentChaosI
{
public:
    TentChaosI() {
        x=1;
        k=1;
    };

    TentChaosI(size_t seed) {
        x=seed&maxVal;
        k=1;
    }

    size_t operator()() {
        k<<=2;
        size_t g=(x+k)&maxVal;
        if((g>>1)<maxVal) {
            x=(g<<1)+1;
        } else {
            x=(maxVal-g)<<1;
        }
        return x;

    }
    static size_t min() {
        return 0;
    }
    static size_t max() {
        return maxVal;
    }
private:
    size_t x;
    size_t k;
    constexpr static const size_t maxVal=~(~size_t(0)<<BitN);

    static_assert (BitN<=64,"BitN mustn't be greater than 64");
    static_assert (BitN>=1,"BitN mustn't be less than 1");
};

}

#endif // TENTCHAOS_H
