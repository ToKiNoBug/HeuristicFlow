// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef Heu_NSGA3_HPP
#define Heu_NSGA3_HPP

#include "NSGA3Base.hpp"

namespace Heu {

/**
 * @brief Parital specialization for NSGA3 using Eigen's array as fitness values
 */
template<typename Var_t,
        size_t ObjNum,
        RecordOption rOpt=DONT_RECORD_FITNESS,
        ReferencePointOption rpOpt=ReferencePointOption::SINGLE_LAYER,
        class Args_t=void,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::initializeFun _iFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::fitnessFun _fFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::crossoverFun _cFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::mutateFun _mFun_=nullptr>
class NSGA3
        : public internal::NSGA3Base<Var_t,ObjNum,rOpt,rpOpt,Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    using Base_t = internal::NSGA3Base<Var_t,ObjNum,rOpt,rpOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_NSGA3ABSTRACT_TYPES

    virtual Fitness_t bestFitness() const {
        Fitness_t best=this->_population.front()._Fitness;
        for(const Gene * i : this->_pfGenes) {
            best=i->_Fitness.min(best);
        }
        return best;
    }

    
    void initializePop() {
        this->makeReferencePoses();
        Base_t::initializePop();
    }
};

}   //  namespace Heu

#endif  //  Heu_NSGA3_HPP
