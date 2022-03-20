// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef Heu_MOGABASE_HPP
#define Heu_MOGABASE_HPP

#include "./MOGAAbstract.hpp"

#ifndef Heu_NO_RTASSERT
#include <assert.h>
#endif

#ifndef Heu_MOGA_MaxRunTimeObjNum
#define Heu_MOGA_MaxRunTimeObjNum 32
#endif

namespace Eigen
{
namespace internal
{
/**
 * @brief MOGA solver base class with compile-time objective count
 * 
 * @tparam Var_t 
 * @tparam ObjNum 
 * @tparam Fitness_t 
 * @tparam fOpt 
 * @tparam rOpt 
 * @tparam Args 
 */
template<typename Var_t,
        int ObjNum,
        FitnessOption fOpt,
        RecordOption rOpt,
        class Args_t,
         typename GAAbstract<Var_t,Eigen::Array<double,ObjNum,1>,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,Eigen::Array<double,ObjNum,1>,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,Eigen::Array<double,ObjNum,1>,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,Eigen::Array<double,ObjNum,1>,Args_t>::mutateFun _mFun_>
class MOGABase
        : public MOGAAbstract<Var_t,ObjNum,fOpt,rOpt,Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
{
    using Base_t = MOGAAbstract<Var_t,ObjNum,fOpt,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
public:
    Heu_MAKE_GABASE_TYPES

    MOGABase() {};
    virtual ~MOGABase() {};

    inline constexpr size_t objectiveNum() const {
        return ObjNum;
    }
};


/**
 * @brief MOGA solver base class with Eigen::Dynamic objective count
 * 
 * @tparam Var_t 
 * @tparam Fitness_t 
 * @tparam fOpt 
 * @tparam rOpt 
 * @tparam Args 
 */
template<typename Var_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        class Args_t,
         typename GAAbstract<Var_t,Eigen::ArrayXd,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,Eigen::ArrayXd,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,Eigen::ArrayXd,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,Eigen::ArrayXd,Args_t>::mutateFun _mFun_>
class MOGABase<Var_t,Eigen::Dynamic,fOpt,rOpt,Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
        : public MOGAAbstract<Var_t,Eigen::Dynamic,fOpt,rOpt,Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
{
public:

    MOGABase() {};
    virtual ~MOGABase() {};

    using Base_t = MOGAAbstract<Var_t,Eigen::Dynamic,fOpt,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_GABASE_TYPES

    inline int objectiveNum() const {
        return _objectiveNum;
    }

    inline void setObjectiveNum(int _objNum) {
#ifndef Heu_NO_RTASSERT
        assert(_objNum>1);
        assert(_objNum<=Heu_MOGA_MaxRunTimeObjNum);
#endif
        _objectiveNum=_objNum;
    }

protected:
    int _objectiveNum;
};

}   //  internal

}   //  namespace Eigen


#endif // MOGABASE_HPP
