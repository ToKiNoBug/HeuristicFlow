

#ifndef HEU_PSOOPTION_HPP
#define HEU_PSOOPTION_HPP

#include <stdint.h>
#include "InternalHeaderCheck.h"
#include <HeuristicFlow/Global>

namespace heu {

// Options to PSO

/**
 * \ingroup HEU_PSO
 * \struct PSOOption
 * \brief Basical option of PSO
 *
 */
struct PSOOption {
 public:
  PSOOption() {
    populationSize = 200;
    maxGeneration = 300;
    maxFailTimes = 100;
    inertiaFactor = 0.8;
    learnFactorP = 2;
    learnFactorG = 2;
  }
  /// size of population, default value is 200
  size_t populationSize;
  /// maximum allowed generation, default value is 300
  size_t maxGeneration;
  /// maximun allowed failtimes, default value is 100
  size_t maxFailTimes;
  /// speed factor, default value is 0.8
  double inertiaFactor;
  /// pBest factor, default value is 2
  double learnFactorP;
  /// gBest factor, default value is 2
  double learnFactorG;
};

}  //  namespace heu

#endif  // HEU_PSOOPTION_HPP
