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

#ifndef HEU_NSGA3_HPP
#define HEU_NSGA3_HPP

#include "InternalHeaderCheck.h"
#include "NSGA3Base.hpp"
#include "DefaultGeneType.hpp"

namespace heu {

/**
 * @brief Parital specialization for NSGA3 using Eigen's array as fitness values
 */

/**
 * \ingroup HEU_GENETIC
 * \brief NSGA3 is the thrid Nondominated Sorting Genetic Algorithm. It's suitable for many
 * objective problems.
 *
 * NSGA3 uses many reference points to maintain a diverse and uniform PF.
 * \sa internal::NSGA3Abstract::select for this special procedure.
 *
 *
 * The default value of template parameters are listed in braces
 * \tparam Var_t Type of decision variable
 * \tparam ObjNum Number of objectives (Eigen::Dynamic for runtime)
 * \tparam rOpt Record fitness or not (don't record)
 * \tparam rpOpt Reference point option (Single layer)
 * \tparam Args_t Other parameters (void)
 * \tparam _iFun_ Initialization function (nullptr)
 * \tparam _fFun_ Fitness function (nullptr)
 * \tparam _cFun_ Crossover function (nullptr)
 * \tparam _mFun_ Mutation function (nullptr)
 *
 * \sa GAOption for ga running parameters
 * \sa SOGA for APIs that all genetic solvers have
 * \sa NSGA2 for APIs that all MOGA solvers have
 *
 * ## APIs that all NSGA3 solvers have:
 * - `const RefMat_t& referencePoints() const` returns a matrix of reference points. Each coloumn is
 * the coordinate of a RP. (Here RP refers to the word reference point)
 * - `size_t referencePointCount() const` number of reference points according to the RP precision.
 *
 *
 * ## APIs that NSGA3 solvers using single-layer RPs have:
 * - `size_t referencePointPrecision() const` returns the RP precision.
 * - `void setReferencePointPrecision(size_t)` set the RP precision
 *
 *
 * ## APIs that NSGA2 solvers using double-layer RPs have:
 * - `size_t innerPrecision() const` returns the precision of inner layer RP.
 * - `size_t outerPrecision() const` returns the precision of outer layer RP.
 * - `void setReferencePointPrecision(size_t i, size_t o)` set the precison of inner and outer
 * layer.
 *
 * \note When using NSGA3 solvers, is strongly recommended to set the precision explicitly before
 * initializing the population. Don't rely on the default value!
 */
template <typename Var_t, int ObjNum, RecordOption rOpt = DONT_RECORD_FITNESS,
          ReferencePointOption rpOpt = ReferencePointOption::SINGLE_LAYER, class Args_t = void,
          typename internal::GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>,
                                        Args_t>::initializeFun _iFun_ = nullptr,
          typename internal::GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::fitnessFun
              _fFun_ = nullptr,
          typename internal::GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>,
                                        Args_t>::crossoverFun _cFun_ = nullptr,
          typename internal::GAAbstract<Var_t, Eigen::Array<double, ObjNum, 1>, Args_t>::mutateFun
              _mFun_ = nullptr>
class NSGA3
    : public internal::NSGA3Base<Var_t, ObjNum, rOpt, rpOpt,
                                 internal::DefaultGene_t<Var_t, Eigen::Array<double, ObjNum, 1>>,
                                 Args_t, _iFun_, _fFun_, _cFun_, _mFun_> {
  using Base_t =
      internal::NSGA3Base<Var_t, ObjNum, rOpt, rpOpt,
                          internal::DefaultGene_t<Var_t, Eigen::Array<double, ObjNum, 1>>, Args_t,
                          _iFun_, _fFun_, _cFun_, _mFun_>;

 public:
  NSGA3() = default;
  ~NSGA3() = default;
  HEU_MAKE_NSGA3ABSTRACT_TYPES(Base_t)

  /**
   * \brief
   *
   */
  inline void initializePop() noexcept {
    this->makeReferencePoses();
    Base_t::initializePop();
  }

  HEU_RELOAD_MEMBERFUCTION_RUN
};

}  //  namespace heu

#endif  //  HEU_NSGA3_HPP
