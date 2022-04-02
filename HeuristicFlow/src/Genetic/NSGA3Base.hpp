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

#ifndef HEU_NSGA3BASE_HPP
#define HEU_NSGA3BASE_HPP

#include "InternalHeaderCheck.h"
#include "NSGA3Abstract.hpp"

namespace heu {

/**
 * \ingroup HEU_GENETIC
 * \brief Different options of reference points.
 *
 * Here's two ways to generate reference points: the Das and Dennis' systematic approach generates a single layer, while
 * Deb and Jain’s method(based on the previous method) generates two layers of RS.
 *
 * The first method places RPs on a M-1 dim hyperplane (where M is the number of objectives) equally inclined on all
 * axes and has a intercept of 1 to all objectives. This method requires relatively higher precision and generates great
 * numbers of RPs for high-dimensional MO problems. This is detrimental to preformance.
 *
 * While the Deb and Jain’s method adds a inner hyperplane(inside layer) whose intercept has smaller value than the
 * boundary layer(normalized hyperplane). Then the algorithm places RPs on two layers with relatively small precision.
 * It works with much fewer RPs.
 */
enum ReferencePointOption {
  SINGLE_LAYER,  ///< Das and Dennis' method
  DOUBLE_LAYER   ///< Deb and Jain's method
};

namespace internal {

/**
 * \ingroup HEU_GENETIC
 * \class NSGA3Base
 * \brief Internal base class for NSGA3(single layer).
 *
 * This class implements different methods to make reference points.
 *
 * \sa NSGA3Base<Var_t, ObjNum, rOpt, DOUBLE_LAYER, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>
 *
 * \tparam Var_t
 * \tparam ObjNum
 * \tparam rOpt
 * \tparam rpOpt Reference point method.
 * \tparam Args_t
 * \tparam _iFun_
 * \tparam _fFun_
 * \tparam _cFun_
 * \tparam _mFun_
 */
template <typename Var_t, int ObjNum, RecordOption rOpt, ReferencePointOption rpOpt, class Args_t,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::initializeFun _iFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun _fFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::crossoverFun _cFun_,
          typename GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun _mFun_>
class NSGA3Base : public NSGA3Abstract<Var_t, ObjNum, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
  using Base_t = NSGA3Abstract<Var_t, ObjNum, rOpt, Args_t, _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  ~NSGA3Base() {}
  HEU_MAKE_NSGA3ABSTRACT_TYPES(Base_t)

  /**
   * \brief Get the precision of RPs
   *
   * \note This function only exists when single layer.
   *
   * \return size_t Precision
   */
  size_t referencePointPrecision() const { return _precision; }

  /**
   * \brief Set the Reference Point Precision object
   *
   * \note This function only exists when single layer.
   *
   * \param p Precison value, not less than 2.
   */
  void setReferencePointPrecision(size_t p) { _precision = p; }

  /**
   * \brief Get the number of RPs
   *
   * \return size_t Number of RPs
   */
  size_t referencePointCount() const { return NchooseK(_precision + this->objectiveNum() - 1, _precision); }

 protected:
  size_t _precision;  ///< Precision of RPs in single layer. Only exists when rpOpt==SINGLE_LAYER

  /**
   * \brief Compute reference points
   *
   * This function is reloaded with different rpOpt.
   */
  void makeReferencePoses() {
    this->referencePoses.resize(this->objectiveNum(), referencePointCount());
    std::vector<Fitness_t> rfP;
    this->computeReferencePointPoses(this->objectiveNum(), _precision, &rfP);
    std::shuffle(rfP.begin(), rfP.end(), global_mt19937());
    for (int c = 0; c < this->referencePoses.cols(); c++) {
      for (int r = 0; r < this->referencePoses.rows(); r++) {
        this->referencePoses(r, c) = rfP[c][r];
      }
    }
  }
};

/**
 * \ingroup HEU_GENETIC
 * \class NSGA3Base
 * \brief Internal base class for NSGA3(double layers).
 *
 * This class implements different methods to make reference points.
 *
 * \sa NSGA3Base
 *
 * \tparam Var_t
 * \tparam ObjNum
 * \tparam rOpt
 * \tparam Args_t
 * \tparam _iFun_
 * \tparam _fFun_
 * \tparam _cFun_
 * \tparam _mFun_
 */
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
  HEU_MAKE_NSGA3ABSTRACT_TYPES(Base_t)

  NSGA3Base() {
    _innerPrecision = 3;
    _outerPrecision = 4;
  };
  ~NSGA3Base() {}

  /**
   * \brief Get the precison of inner layer
   *
   * \note This function only exists with double layers.
   *
   * \return size_t Precision of inner layer
   */
  size_t innerPrecision() const { return _innerPrecision; }

  /**
   * \brief Get the precison of outer layer
   *
   * \note This function only exists with double layers.
   *
   * \return size_t Precision of outer layer
   */
  size_t outerPrecision() const { return _outerPrecision; }

  /**
   * \brief Set the precison of inner and outer layers
   *
   * \note This function only exists with double layers.
   *
   * \param i Inner layer precision
   * \param o Outer layer precision
   */
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
  size_t _innerPrecision;  ///< Precision of inner layer
  size_t _outerPrecision;  ///< Precision of outer layer

  void makeReferencePoses() {
    this->referencePoses.resize(this->objectiveNum(), referencePointCount());
    std::vector<Fitness_t> irfP, orfP;
    std::shuffle(irfP.begin(), irfP.end(), global_mt19937());
    std::shuffle(orfP.begin(), orfP.end(), global_mt19937());
    this->computeReferencePointPoses(this->objectiveNum(), _innerPrecision, &irfP);
    this->computeReferencePointPoses(this->objectiveNum(), _outerPrecision, &orfP);

    for (size_t c = 0; c < (size_t)this->referencePoses.cols(); c++) {
      for (size_t r = 0; r < this->objectiveNum(); r++) {
        if (c < irfP.size()) {
          this->referencePoses(r, c) = irfP[c][r] * 0.70710678118654752440;
        } else {
          this->referencePoses(r, c) = orfP[c - irfP.size()][r];
        }
      }
    }

  }  //  makeReferencePoses()
};

}  //  namespace internal

}  //  namespace heu

#endif  //  HEU_NSGA3BASE_HPP
