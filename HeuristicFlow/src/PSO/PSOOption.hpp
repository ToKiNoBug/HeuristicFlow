// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_PSOOPTION_HPP
#define EIGEN_HEU_PSOOPTION_HPP

#include <stdint.h>
#include "InternalHeaderCheck.h"
#include "../../Global"

/**
 * \defgroup HEU_PSO Particle Swarm Optimization algorithm
 *
 * This module provides fundamental implementation of PSO.
 */

namespace Eigen {

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

}  //  namespace Eigen

#endif  // EIGEN_HEU_PSOOPTION_HPP
