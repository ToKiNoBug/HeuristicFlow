// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_NSGA3BASE_HPP
#define EIGEN_HEU_NSGA3BASE_HPP

#include "NSGA3Abstract.hpp"

namespace Eigen {

/**
 * @brief Layers of Reference points
 *
 */
enum ReferencePointOption { SINGLE_LAYER, DOUBLE_LAYER };

namespace internal {

template <typename Var_t, int ObjNum, RecordOption rOpt, ReferencePointOption rpOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_>
class NSGA3Base : public NSGA3Abstract<Var_t, ObjNum, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
  using Base_t = NSGA3Abstract<Var_t, ObjNum, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  EIGEN_HEU_MAKE_NSGA3ABSTRACT_TYPES(Base_t)

  size_t referencePointPrecision() const { return _precision; }

  void setReferencePointPrecision(size_t p) { _precision = p; }

  size_t referencePointCount() const { return NchooseK(_precision + this->objectiveNum() - 1, _precision); }

 protected:
  size_t _precision;

  void makeReferencePoses() {
    this->referencePoses.resize(this->objectiveNum(), referencePointCount());
    std::vector<Fitness_t> rfP;
    this->computeReferencePointPoses(this->objectiveNum(), _precision, &rfP);
    std::shuffle(rfP.begin(), rfP.end(), global_mt19937);
    for (size_t c = 0; c < this->referencePoses.cols(); c++) {
      for (size_t r = 0; r < this->referencePoses.rows(); r++) {
        this->referencePoses(r, c) = rfP[c][r];
      }
    }
  }
};

template <typename Var_t, int ObjNum, RecordOption rOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_>
class NSGA3Base<Var_t, ObjNum, rOpt, DOUBLE_LAYER, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>
    : public NSGA3Abstract<Var_t, ObjNum, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
 private:
  using Base_t = NSGA3Abstract<Var_t, ObjNum, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  EIGEN_HEU_MAKE_NSGA3ABSTRACT_TYPES(Base_t)

  NSGA3Base() {
    _innerPrecision = 3;
    _outerPrecision = 4;
  };

  size_t innerPrecision() const { return _innerPrecision; }

  size_t outerPrecision() const { return _outerPrecision; }

  void setReferencePointPrecision(size_t i, size_t o) {
    assert(i >= 2);
    assert(o >= 2);
    _innerPrecision = i;
    _outerPrecision = o;
  }

  size_t referencePointCount() const {
    return NchooseK(_innerPrecision + this->objectiveNum() - 1, _innerPrecision) +
           NchooseK(_outerPrecision + this->objectiveNum() - 1, _outerPrecision);
  }

 protected:
  size_t _innerPrecision;
  size_t _outerPrecision;

  void makeReferencePoses() {
    this->referencePoses.resize(this->objectiveNum(), referencePointCount());
    std::vector<Fitness_t> irfP, orfP;
    std::shuffle(irfP.begin(), irfP.end(), global_mt19937());
    std::shuffle(orfP.begin(), orfP.end(), global_mt19937());
    this->computeReferencePointPoses(this->objectiveNum(), _innerPrecision, &irfP);
    this->computeReferencePointPoses(this->objectiveNum(), _outerPrecision, &orfP);

    for (int c = 0; c < this->referencePoses.cols(); c++) {
      for (int r = 0; r < this->objectiveNum(); r++) {
        if (c < irfP.size()) {
          this->referencePoses(r, c) = irfP[c][r] * M_SQRT1_2;
        } else {
          this->referencePoses(r, c) = orfP[c - irfP.size()][r];
        }
      }
    }

  }  //  makeReferencePoses()
};

}  //  namespace internal

}  //  namespace Eigen

#endif  //  EIGEN_HEU_NSGA3BASE_HPP
