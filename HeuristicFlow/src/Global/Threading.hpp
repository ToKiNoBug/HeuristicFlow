/*
 Copyright © 2021-2022  TokiNoBug
This file is part of HeuristicFlow.

    HeuristicFlow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HeuristicFlow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HeuristicFlow.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef HEU_THREADING_HPP
#define HEU_THREADING_HPP

#ifdef EIGEN_HAS_OPENMP
#include <omp.h>
#endif  //  EIGEN_HAS_OPENMP

#include <thread>

#include <assert.h>

#include "InternalHeaderCheck.h"

#include "Macros.hpp"

namespace heu {

int threadNum() noexcept;
inline void setThreadNum(int _thN) noexcept;

namespace {

class thNumWrapper {
  thNumWrapper() = delete;
  friend int ::heu::threadNum() noexcept;
  friend void ::heu::setThreadNum(int) noexcept;
  static inline int& internalThreadNum() noexcept {
    static int thNum = std::thread::hardware_concurrency();
    return thNum;
  }
};

}  // namespace

/**
 * \ingroup HEU_GLOBAL
 * \brief Query the number of threads that will be used. The default value euqals to
 * `std::thread::hardware_concurrency()`.
 *
 * \return int Number of threads that will be used.
 */
inline int threadNum() noexcept { return thNumWrapper::internalThreadNum(); }

/**
 * \ingroup HEU_GLOBAL
 * \brief Set the number of threads that will be used. It will cause a assertion failure if `_thN`
 * is not a positive integer.
 *
 * \param _thN Number of threads that will be used.
 */
inline void setThreadNum(int _thN) noexcept {
  HEU_ASSERT(_thN > 0);
  thNumWrapper::internalThreadNum() = _thN;
}

// if HEU_NO_THREADS is defined, HEU_HAS_OPENMP won't be defined, preventing following algorithms to
// apply multi-threading. Otherwise if Eigen has OpenMP, then HeuristicFlow will use it as well.
#ifndef HEU_NO_THREADS

#ifdef EIGEN_HAS_OPENMP

#ifndef HEU_HAS_OPENMP
#define HEU_HAS_OPENMP
#endif  //  HEU_HAS_OPENMP

#endif  //  EIGEN_HAS_OPENMP

#endif  //  HEU_NO_THREADS

}  // namespace heu

#endif  // HEU_THREADING_HPP
