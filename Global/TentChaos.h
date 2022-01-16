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
