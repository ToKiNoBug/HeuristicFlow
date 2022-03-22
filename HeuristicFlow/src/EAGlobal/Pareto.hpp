// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_PARETO_HPP
#define EIGEN_HEU_PARETO_HPP

#include "InternalHeaderCheck.h"
#include "../../Global"

/**
 * \defgroup HEU_EAGlobal EAGlobal
 *
 */

namespace Eigen {

namespace internal {

/**
 * \ingroup HEU_EAGlobal
 * \struct Pareto
 * \brief Pareto optimality in multi-objective problems
 *
 * \tparam ObjNum Number of objectives, Eigen::Dynamic for runtime objs
 * \tparam fOpt Whether greater or less fitness means better.
 */
template <int ObjNum, FitnessOption fOpt>
struct Pareto {
  static_assert(ObjNum > 0 || ObjNum == Eigen::Dynamic, "ObjNum should be positive or dynamic(-1)");
  static_assert(ObjNum != 1, "You assigned 1 objective for multi-objective problem");

  using Fitness_t = Eigen::Array<double, ObjNum, 1>;
  /**
   * \brief Whether A dominates B
   *
   * \param A Pointer to the first fitness value
   * \param B Pointer to the second fitness value
   * \return true A is not worse than B on all objectives AND A has at least one objective that is better than B.
   * \return false A doesn't dominate B
   */
  static bool isStrongDominate(const Fitness_t* A, const Fitness_t* B) {
    bool isNotWorse, isBetter;
    if (fOpt == FITNESS_GREATER_BETTER) {
      isNotWorse = ((*A) >= (*B)).all();
      isBetter = ((*A) > (*B)).any();
    } else {
      isNotWorse = ((*A) <= (*B)).all();
      isBetter = ((*A) < (*B)).any();
    }

    return isNotWorse && isBetter;
  }
};

}  //    namespace internal
}  //    namespace Eigen

#endif  // EIGEN_HEU_PARETO_HPP
