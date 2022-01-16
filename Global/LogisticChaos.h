#ifndef LOGISTICCHAOS_H
#define LOGISTICCHAOS_H

#include <cstdint>

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

    void iterate(size_t It) {
        for(size_t i=0;i<It;i++) {
            val*=miu*(1-val);
        }
    }
private:
    double val;
    constexpr static const double miu=4.0;
};

#endif // LOGISTICCHAOS_H
