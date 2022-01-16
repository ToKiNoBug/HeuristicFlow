#ifndef LOGISTICCHAOS_H
#define LOGISTICCHAOS_H

#include <cstdint>
namespace OptimT {

class LogisticChaos
{
public:
    LogisticChaos() {
        val=0.1;
    }

    LogisticChaos(double seed) {
        val=seed;
    }

    double operator()() {
        val*=miu*(1-val);
        return val;
    }

    double value() const {
        return val;
    }

    void makeSequence(double *dst,size_t L) {
        if(L<=0) return;
        dst[0]=operator()();
        for(size_t i=1;i<L;i++) {
            dst[i]=miu*dst[i-1]*(1-dst[i-1]);
        }
        val=dst[L-1];
    }

    void iterate(size_t It) {
        for(size_t i=0;i<It;i++) {
            val*=miu*(1-val);
        }
    }
private:
    double val;
    constexpr static const double miu=4.0;
};

}
#endif // LOGISTICCHAOS_H
