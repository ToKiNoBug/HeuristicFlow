/*
 Copyright © 2022  TokiNoBug
This file is part of Heuristic.

    Heuristic is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Heuristic is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Heuristic.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef Heu_PSOPARAMETERPACK_HPP
#define Heu_PSOPARAMETERPACK_HPP

#include "../../Global"

namespace Heu
{


namespace HeuPrivate
{

Heu_MAKE_FUNAREA(iFun,iFun,PSO)
Heu_MAKE_FUNAREA(fFun,fFun,PSO)

}

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
        HeuPrivate::iFunArea_PSO<void,Var_t *,Var_t *,
            const Var_t *,const Var_t *,
            const Var_t *,const Args_t *>::template funBody<_i>;

    template <fFun_t _i>
    using fFunBody = typename
        HeuPrivate::fFunArea_PSO<void,const Var_t * ,
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
        HeuPrivate::iFunArea_PSO<void,Var_t *,Var_t *,
            const Var_t *,const Var_t *,
            const Var_t *>::template funBody<_i>;

    template <fFun_t _i>
    using fFunBody = typename
        HeuPrivate::fFunArea_PSO<void,const Var_t *,Fitness_t *>::template funBody<_i>;


    static const bool HasParameters=false;
};

#define Heu_MAKE_PSOPARAMETERPACK_TYPES \
using Args_t = typename Base_t::Args_t; \
using iFun_t = typename Base_t::iFun_t; \
using fFun_t = typename Base_t::fFun_t;

} // namespace Heu


#endif  //  Heu_PSOPARAMETERPACK_HPP
