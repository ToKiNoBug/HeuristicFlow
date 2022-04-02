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

#ifndef HEU_PSODEFAULTS_HPP
#define HEU_PSODEFAULTS_HPP

#include "InternalHeaderCheck.h"
#include "PSO.hpp"

namespace heu {

namespace internal {
template <typename Var_t, bool isVarEigenTypes>
struct __impl_PSODefaults {
  // Candidate function for initialization
  inline static void impl_iFun(Var_t *x, Var_t *v, const Var_t *xMin, const Var_t *xMax, const Var_t *) {
    for (int idx = 0; idx < xMin->size(); idx++) {
      x->operator[](idx) = ei_randD(xMin->operator[](idx), xMax->operator[](idx));
      v->operator[](idx) = 0;
    }
  }
};

template <typename Var_t>
struct __impl_PSODefaults<Var_t, true> {
  // Candidate function for initialization
  inline static void impl_iFun(Var_t *x, Var_t *v, const Var_t *xMin, const Var_t *xMax, const Var_t *) {
    x->setRandom(xMin->size(), 1);
    (*x) *= (*xMax - *xMin) / 2;
    (*x) += (*xMin + *xMax) / 2;
    v->setZero(xMin->size(), 1);
  }
};
}  // namespace internal

/**
 * \ingroup HEU_PSO
 * \struct PSODefaults
 * \brief PSODefaults provides candidate functions for initialization
 *
 * This struct has a specialization for `Args_t` is `void`.
 *
 * \sa PSODefaults<Var_t, isVarEigenTypes, void>
 *
 * \tparam Var_t Type of decision variable
 * \tparam isVarEigenTypes If `Var_t` is Eigen's Array/Matrix(s).
 * \tparam Args_t Pseudo-global variables type in solver.
 */
template <typename Var_t, bool isVarEigenTypes = false, typename Args_t = void>
struct PSODefaults {
  static_assert(!std::is_same<Args_t, void>::value, "Wrong specialization of PSODefaults");

  /**
   * \brief Default initialize function with args
   *
   * \param pos Position to be initialized
   * \param velocity Velocity to be initialized
   * \param pMin Minmum position
   * \param pMax Maximum position
   * \param vMax Maximum speed (absolute value)
   *
   * \sa PSODefaults<Var_t, isVarEigenTypes, void>::iFun
   */
  inline static void iFun(Var_t *pos, Var_t *velocity, const Var_t *pMin, const Var_t *pMax, const Var_t *vMax,
                          const Args_t *) {
    internal::__impl_PSODefaults<Var_t, isVarEigenTypes>::impl_iFun(pos, velocity, pMin, pMax, vMax);
  }
};

/**
 * \ingroup HEU_PSO
 * \struct PSODefaults<Var_t, isVarEigenTypes, void>
 * \brief Specialization when `Args_t` is `void`.
 *
 * \sa PSODefaults
 *
 * \tparam Var_t
 * \tparam isVarEigenTypes
 */
template <typename Var_t, bool isVarEigenTypes>
struct PSODefaults<Var_t, isVarEigenTypes, void> {
  /**
   * \brief Default initialize function without args
   *
   * \param pos Position to be initialized
   * \param velocity Velocity to be initialized
   * \param pMin Minmum position
   * \param pMax Maximum position
   * \param vMax Maximum speed (absolute value)
   *
   * \sa PSODefaults::iFun
   */
  inline static void iFun(Var_t *pos, Var_t *velocity, const Var_t *pMin, const Var_t *pMax, const Var_t *vMax) {
    internal::__impl_PSODefaults<Var_t, isVarEigenTypes>::impl_iFun(pos, velocity, pMin, pMax, vMax);
  }
};

}  // namespace heu

#endif  //  HEU_PSODEFAULTS_HPP
