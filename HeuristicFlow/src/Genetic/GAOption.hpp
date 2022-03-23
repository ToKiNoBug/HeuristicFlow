// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_GAOPTION_HPP
#define EIGEN_HEU_GAOPTION_HPP

#include <stdint.h>

#include "InternalHeaderCheck.h"

namespace Eigen {

/// options about GA algorithm

/**
 * \ingroup HEU_Genetic
 * \struct GAOption
 * \brief GAOption is a struct encapsulating several general options for most types of genetic algorithm.
 *
 * GAOption is a non-template struct, which means that it's suitable to most genetic algorithm. It's a simple struct
 * with all members public and without any function except constructor.
 */
struct GAOption {
 public:
  /**
   * \brief Construct and initialize all members to their default values.
   *
   */
  GAOption() {
    populationSize = 100;
    maxFailTimes = 50;
    maxGenerations = 300;
    crossoverProb = 0.8;
    mutateProb = 0.05;
  }

  /**
   * \brief Size of the population. Default value is 100
   *
   */
  size_t populationSize;

  /**
   * \brief GA will stop once best solution hasn't been improved for continuous maxFailTimes generations. Default value
   * is 100.
   *
   */
  size_t maxFailTimes;

  /**
   * \brief Maximum generation. GA will stop once reached this limitation. Default value is 300
   *
   */
  size_t maxGenerations;

  /**
   * \brief Probability of a non-elite individual to join crossover. Default value is 80%
   *
   */
  double crossoverProb;

  /**
   * \brief Probability of a non-elite individual to get mutated. Default value is 5%
   *
   */
  double mutateProb;
};

}  //    namespace Eigen

#endif  // EIGEN_HEU_GAOPTION_HPP
