// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_GAABSTRACT_HPP
#define EIGEN_HEU_GAABSTRACT_HPP

#include <type_traits>

#include "../../Global"

namespace Eigen
{

namespace internal
{

EIGEN_HEU_MAKE_FUNAREA(iFun,iFun,GA)
EIGEN_HEU_MAKE_FUNAREA(fFun,fFun,GA)
EIGEN_HEU_MAKE_FUNAREA(cFun,cFun,GA)
EIGEN_HEU_MAKE_FUNAREA(mFun,mFun,GA)


/**
 * @brief The GAAbstract class declares 4 operator functions and related reloaded functions.
 */
template<typename Var_t,typename Fitness_t,class Args_t>
class GAAbstract
{
public:
   GAAbstract() {};
    virtual ~GAAbstract() {};

   ///Function to initialize Var
   using initializeFun = void(*)(Var_t*,const Args_t*);
   ///Function to calculate fitness for Var
   using fitnessFun = void(*)(const Var_t*,const Args_t*,Fitness_t*);
   ///Function to apply crossover for Var
   using crossoverFun = void(*)(const Var_t*,const Var_t*,Var_t*,Var_t*,const Args_t*);
   ///Function to apply mutate for Var
   using mutateFun = void(*)(const Var_t *,Var_t * ,const Args_t *);

   using ArgsType = Args_t;

   template<initializeFun i>
   using iFunBody = typename iFunArea_GA<void,Var_t*,const Args_t *>::template funBody<i>;

   template<fitnessFun f>
   using fFunBody = typename fFunArea_GA<void,const Var_t*,const Args_t *,Fitness_t*>::template funBody<f>;

   template<crossoverFun c>
   using cFunBody = typename cFunArea_GA<void,const Var_t*,const Var_t*,Var_t*,Var_t*,const Args_t *>::template funBody<c>;

   template<mutateFun m>
   using mFunBody= typename mFunArea_GA<void,const Var_t *,Var_t *,const Args_t *>::template funBody<m>;

    const Args_t & args() const {
        return _args;
    }

    void setArgs(const Args_t & a) {
        _args=a;
    }

    static const bool HasParameters=true;

protected:
Args_t _args;
};


template<typename Var_t,typename Fitness_t>
class GAAbstract<Var_t,Fitness_t,void>
{
public:
    GAAbstract() {};
    virtual ~GAAbstract() {};

    ///Function to initialize Var
    using initializeFun = void(*)(Var_t*);
    ///Function to calculate fitness for Var
    using fitnessFun = void(*)(const Var_t*,Fitness_t*);
    ///Function to apply crossover for Var
    using crossoverFun = void(*)(const Var_t*,const Var_t*,Var_t*,Var_t*);
    ///Function to apply mutate for Var
    using mutateFun = void(*)(const Var_t * ,Var_t*);

    using ArgsType = void;

    template<initializeFun i>
    using iFunBody = typename iFunArea_GA<void,Var_t*>::template funBody<i>;

    template<fitnessFun f>
    using fFunBody = typename fFunArea_GA<void,const Var_t*,Fitness_t*>::template funBody<f>;

    template<crossoverFun c>
    using cFunBody = typename cFunArea_GA<void,const Var_t*,const Var_t*,Var_t*,Var_t*>::template funBody<c>;

    template<mutateFun m>
    using mFunBody= typename mFunArea_GA<void,const Var_t *,Var_t *>::template funBody<m>;

    static const bool HasParameters=false;
};

#define EIGEN_HEU_MAKE_GAABSTRACT_TYPES(Base_t) \
using initializeFun = typename Base_t::initializeFun; \
using fitnessFun = typename Base_t::fitnessFun; \
using crossoverFun = typename Base_t::crossoverFun; \
using mutateFun = typename Base_t::mutateFun; \
using ArgsType = typename Base_t::ArgsType;

}   //  internal

}   //  namespace Eigen

#endif  //  EIGEN_HEU_GAABSTRACT_HPP
