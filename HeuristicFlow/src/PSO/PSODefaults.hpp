// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef EIGEN_HEU_PSODEFAULTS_HPP
#define EIGEN_HEU_PSODEFAULTS_HPP
#include "PSO.hpp"

namespace Eigen
{

namespace internal
{
    template<typename Var_t,bool isVarEigenTypes>
    struct __impl_PSODefaults
    {
    ///Candidate function for initialization
    inline static void impl_iFun(Var_t *x,Var_t *v,
    const Var_t * xMin,const Var_t * xMax,
    const Var_t *) {
        for(int idx=0;idx<xMin->size();idx++) {
            x->operator[](idx)=randD(xMin->operator[](idx),xMax->operator[](idx));
            v->operator[](idx)=0;
        }
    }
    };


    template<typename Var_t>
    struct __impl_PSODefaults<Var_t,true>
    {

    ///Candidate function for initialization
    inline static void impl_iFun(Var_t *x,Var_t *v,
        const Var_t * xMin,const Var_t * xMax,
        const Var_t *) {
        x->setRandom(xMin->size(),1);
        (*x)*=(*xMax-*xMin)/2;
        (*x)+=(*xMin+*xMax)/2;
        v->setZero(xMin->size(),1);
    }
    };
}


template<typename Var_t,bool isVarEigenTypes=false,typename Args_t=void>
struct PSODefaults
{
    static_assert(!std::is_same<Args_t,void>::value,"Wrong specialization of PSODefaults");

    inline static void iFun(Var_t * pos,Var_t * velocity,
        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax,const Args_t *) 
    {
        internal::__impl_PSODefaults<Var_t,isVarEigenTypes>::impl_iFun(pos,velocity,pMin,pMax,vMax);
    }
};


template<typename Var_t,bool isVarEigenTypes>
struct PSODefaults<Var_t,isVarEigenTypes,void>
{
    inline static void iFun(Var_t * pos,Var_t * velocity,
        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax) 
    {
        internal::__impl_PSODefaults<Var_t,isVarEigenTypes>::impl_iFun(pos,velocity,pMin,pMax,vMax);
    }
};


}   //  namespace Eigen

#endif  //  EIGEN_HEU_PSODEFAULTS_HPP
