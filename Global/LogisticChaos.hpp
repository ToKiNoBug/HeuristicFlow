#ifndef LOGISTICCHAOS_H
#define LOGISTICCHAOS_H

#include <cstdint>
#include "Macros.hpp"

namespace OptimT {

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
