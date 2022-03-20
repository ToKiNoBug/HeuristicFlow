// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_PSOBASE_HPP
#define Heu_PSOBASE_HPP

#include "PSOOption.hpp"
#include "PSOAbstrcat.hpp"

namespace Eigen
{

namespace internal
{

///Abstrcat base class for most PSO solvers
template<class Var_t,int DIM,class Fitness_t,RecordOption Record,class Arg_t,
         typename PSOParameterPack<Var_t,Fitness_t,Arg_t>::iFun_t _iFun_,
         typename PSOParameterPack<Var_t,Fitness_t,Arg_t>::fFun_t _fFun_>
class PSOBase : public PSOAbstract<Var_t,Fitness_t,Record,Arg_t,_iFun_,_fFun_>
{
    using Base_t = PSOAbstract<Var_t,Fitness_t,Record,Arg_t,_iFun_,_fFun_>;
public:
    Heu_MAKE_PSOABSTRACT_TYPES

public:
    PSOBase() {};
    virtual ~PSOBase() {};

    static constexpr int dimensions() {
        return DIM;
    }

protected:
    static const int dims=DIM;

};

///partial specialization for PSOBase with Runtime dimensions
template<class Var_t,class Fitness_t,RecordOption Record,class Arg_t,
         typename PSOParameterPack<Var_t,Fitness_t,Arg_t>::iFun_t _iFun_,
         typename PSOParameterPack<Var_t,Fitness_t,Arg_t>::fFun_t _fFun_>
class PSOBase<Var_t,Eigen::Dynamic,Fitness_t,Record,Arg_t,_iFun_,_fFun_>
        : public PSOAbstract<Var_t,Fitness_t,Record,Arg_t,_iFun_,_fFun_>
{
    using Base_t = PSOAbstract<Var_t,Fitness_t,Record,Arg_t,_iFun_,_fFun_>;
public:
    Heu_MAKE_PSOABSTRACT_TYPES

    inline int dimensions() const {
        return dims;
    }

    inline void setDimensions(int d) {
        dims=d;
    }

protected:
    int dims;

};

}   //  internal

}   //  namespace Eigen
#endif // PSOBASE_HPP
