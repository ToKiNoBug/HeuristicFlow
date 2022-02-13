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

#ifndef Heu_LOGISTICCHAOS_HPP
#define Heu_LOGISTICCHAOS_HPP

#include <cstdint>
#include "Macros.hpp"

namespace Heu {

class LogisticChaos
{
public:
    LogisticChaos() {
        _value=0.1;
    }

    LogisticChaos(double seed) {
        _value=seed;
    }

    double operator()() {
        _value*=miu*(1-_value);
        return _value;
    }

    double value() const {
        return _value;
    }

    void makeSequence(double *dst,size_t L) {
        if(L<=0) return;
        dst[0]=operator()();
        for(size_t i=1;i<L;i++) {
            dst[i]=miu*dst[i-1]*(1-dst[i-1]);
        }
        _value=dst[L-1];
    }

    void iterate(size_t It) {
        for(size_t i=0;i<It;i++) {
            _value*=miu*(1-_value);
        }
    }
private:
    double _value;
    constexpr static const double miu=4.0;
};

}
#endif // LOGISTICCHAOS_H
