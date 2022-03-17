// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_TENTCHAOS_H
#define Heu_TENTCHAOS_H

#include <stdint.h>

namespace Heu {

template <uint8_t BitN>
class TentChaosI
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
