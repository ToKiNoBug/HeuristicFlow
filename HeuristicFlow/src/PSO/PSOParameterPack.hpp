// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef EIGEN_HEU_PSOPARAMETERPACK_HPP
#define EIGEN_HEU_PSOPARAMETERPACK_HPP

#include "../../Global"

namespace Eigen
{


namespace internal
{

Heu_MAKE_FUNAREA(iFun,iFun,PSO)
Heu_MAKE_FUNAREA(fFun,fFun,PSO)


template<class Var_t,class Fitness_t,class Arg_t=void>
class PSOParameterPack
{
public:
    PSOParameterPack() {};
    virtual ~PSOParameterPack() {};
    using Args_t = Arg_t;

    using iFun_t = void(*)(Var_t * pos,Var_t * velocity,
        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax,const Args_t *);
    using fFun_t = void(*)(const Var_t * ,const Args_t *,Fitness_t *);


    template <iFun_t _i>
    using iFunBody = typename
        iFunArea_PSO<void,Var_t *,Var_t *,
            const Var_t *,const Var_t *,
            const Var_t *,const Args_t *>::template funBody<_i>;

    template <fFun_t _i>
    using fFunBody = typename
        fFunArea_PSO<void,const Var_t * ,
            const Args_t *,Fitness_t *>::template funBody<_i>;


    void setArgs(const Arg_t & a) {
        _arg=a;
    }

    const Arg_t & args() const {
        return _arg;
    }

    static const bool HasParameters=true;

protected:
    Arg_t _arg;
};


template<class Var_t,class Fitness_t>
class PSOParameterPack<Var_t,Fitness_t,void>
{
public:
    PSOParameterPack() {};
    virtual ~PSOParameterPack() {};
    using Args_t = void;

    using iFun_t = void(*)(Var_t * pos,Var_t * velocity,
        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax);
    using fFun_t = void(*)(const Var_t * ,Fitness_t *);

    template <iFun_t _i>
    using iFunBody = typename
        iFunArea_PSO<void,Var_t *,Var_t *,
            const Var_t *,const Var_t *,
            const Var_t *>::template funBody<_i>;

    template <fFun_t _i>
    using fFunBody = typename
        fFunArea_PSO<void,const Var_t *,Fitness_t *>::template funBody<_i>;


    static const bool HasParameters=false;
};

#define Heu_MAKE_PSOPARAMETERPACK_TYPES \
using Args_t = typename Base_t::Args_t; \
using iFun_t = typename Base_t::iFun_t; \
using fFun_t = typename Base_t::fFun_t;

}   //  internal

} // namespace Eigen


#endif  //  EIGEN_HEU_PSOPARAMETERPACK_HPP
