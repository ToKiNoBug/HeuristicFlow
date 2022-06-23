#ifndef HEU_AOSOPTION_HPP
#define HEU_AOSOPTION_HPP

#include <stdint.h>

#include "InternalHeaderCheck.hpp"

namespace heu {

struct AOSOption {
  AOSOption() {
    electronNum = 50;
    maxGeneration = 30;
    maxEarlyStop = 20;
    photonRate = 0.1;
  }

  ~AOSOption() = default;

  size_t electronNum;
  size_t maxGeneration;
  size_t maxEarlyStop;
  double photonRate;
};

}  // namespace heu

#endif  //  HEU_AOSOPTION_HPP