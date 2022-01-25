#ifndef CHAOTIC_HPP
#define CHAOTIC_HPP
#include "./LogisticChaos.hpp"
namespace OptimT {

extern LogisticChaos global_logistic;

///logistic random number in range (0,1)
inline double logisticD() {
    return global_logistic();
}

}   // OptimT

#endif // CHAOTIC_HPP