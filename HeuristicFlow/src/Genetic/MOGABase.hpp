// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_MOGABASE_HPP
#define EIGEN_HEU_MOGABASE_HPP

#include <assert.h>

#include "InternalHeaderCheck.h"
#include "MOGAAbstract.hpp"

#ifndef Heu_MOGA_MaxRunTimeObjNum
#define Heu_MOGA_MaxRunTimeObjNum 32
#endif

namespace Eigen {
namespace internal {

/**
 * \ingroup HEU_Genetic
 * \class MOGABase
 * \brief Base class for multiple objective genetic solvers.
 *
 * MOGABase maintains the numbers of objectives.
 *
 * This step of inheritance aims to support solvers with fixed and dynamic objective numbers. This template class has a
 * default implementation for fixed objective numbers and a partial specialization for dynamic. Lots of code-copying can
 * be avoided by such inheritance and parital specialization.
 *
 * \tparam Var_t
 * \tparam ObjNum
 * \tparam fOpt
 * \tparam rOpt
 * \tparam Args_t
 * \tparam _iFun_
 * \tparam _fFun_
 * \tparam _cFun_
 * \tparam _mFun_
 */
template <typename Var_t, int ObjNum, FitnessOption fOpt, RecordOption rOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_>
class MOGABase : public MOGAAbstract<Var_t, ObjNum, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
  using Base_t = MOGAAbstract<Var_t, ObjNum, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  EIGEN_HEU_MAKE_GABASE_TYPES(Base_t)

  /**
   * \brief Returns numbers of objectives
   *
   * \return constexpr size_t Numbers of objecvites at compile time.
   */
  inline constexpr size_t objectiveNum() const { return ObjNum; }
};

/**
 * \ingroup HEU_Genetic
 * \class MOGABase
 * \brief Base class for multiple objective genetic solvers.
 *
 * \tparam Var_t
 * \tparam fOpt
 * \tparam rOpt
 * \tparam Args_t
 * \tparam _iFun_
 * \tparam _fFun_
 * \tparam _cFun_
 * \tparam _mFun_
 */
template <typename Var_t, FitnessOption fOpt, RecordOption rOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::ArrayXd, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::ArrayXd, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::ArrayXd, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::ArrayXd, Args_t>::mutateFun _mFun_>
class MOGABase<Var_t, Eigen::Dynamic, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>
    : public MOGAAbstract<Var_t, Eigen::Dynamic, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
 public:
  using Base_t = MOGAAbstract<Var_t, Eigen::Dynamic, fOpt, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;
  EIGEN_HEU_MAKE_GABASE_TYPES(Base_t)

  MOGABase() { _objectiveNum = 0; }

  /**
   * \brief Return the value of objecvites
   *
   * \return int Number of objectives.
   */
  inline int objectiveNum() const { return _objectiveNum; }

  /**
   * \brief Set the Objective Num object
   *
   * \note Runtime assertion will fail if given _objNum is less than 2, or it exceeds Heu_MOGA_MaxRunTimeObjNum.
   * \note This member function exists only when template parameter ObjNum is Eigen::Dynamic.
   *
   * \param _objNum Number of objectives
   */
  inline void setObjectiveNum(int _objNum) {
#ifndef Heu_NO_RTASSERT
    assert(_objNum > 1);
    assert(_objNum <= Heu_MOGA_MaxRunTimeObjNum);
#endif
    _objectiveNum = _objNum;
  }

 protected:
  int _objectiveNum;  ///< Default value is 0
};

}  //  namespace internal

}  //  namespace Eigen

#endif  // EIGEN_HEU_ MOGABASE_HPP
