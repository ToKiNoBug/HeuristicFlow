// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef Heu_NSGA2_HPP
#define Heu_NSGA2_HPP

#include "./NSGA2Base.hpp"
#include <type_traits>
namespace Heu
{


/**
 * @brief NSGA2 MOGA solver. Suitable for not too many objectives.
 * 
 * @tparam Var_t Type of decisition variable.
 * @tparam ObjNum Numbers of objectives.
 * @tparam isGreaterBetter Whether greater fitness value means better.
 * @tparam Record Whether the solver records fitness changelog.
 * @tparam Args_t Type of other parameters.
 * @tparam _iFun_ Compile-time iFun, use nullptr for runtime
 * @tparam _fFun_ Compile-time fFun, use nullptr for runtime
 * @tparam _cFun_ Compile-time cFun, use nullptr for runtime
 * @tparam _mFun_ Compile-time mFun, use nullptr for runtime
 */
template<typename Var_t,
         size_t ObjNum,
         FitnessOption isGreaterBetter=FITNESS_LESS_BETTER,
         RecordOption Record=DONT_RECORD_FITNESS,
         class Args_t=void,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::initializeFun _iFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::fitnessFun _fFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::crossoverFun _cFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::mutateFun _mFun_=nullptr>
class NSGA2
    : public internal::NSGA2Base<Var_t,
                    ObjNum,
                    isGreaterBetter,
                    Record,
                    Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    using Base_t = internal::NSGA2Base<Var_t,
                    ObjNum,
                    isGreaterBetter,
                    Record,
                    Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_NSGABASE_TYPES

    using infoUnit2 = typename Base_t::infoUnit2;

    NSGA2() {}

    virtual ~NSGA2() {}
    
protected:

private:

};

}   //  Heu

#endif // NSGA2_HPP
