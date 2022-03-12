/*
 Copyright © 2022  TokiNoBug
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

#ifndef Heu_GAABSTRACT_HPP
#define Heu_GAABSTRACT_HPP

#include <type_traits>

#include "../../Global"

namespace Heu {

namespace HeuPrivate {

Heu_MAKE_FUNAREA(iFun,iFun,GA)
Heu_MAKE_FUNAREA(fFun,fFun,GA)
Heu_MAKE_FUNAREA(cFun,cFun,GA)
Heu_MAKE_FUNAREA(mFun,mFun,GA)

}

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
   using mutateFun = initializeFun;

   using ArgsType = Args_t;

   template<initializeFun i>
   using iFunBody = typename HeuPrivate::iFunArea_GA<void,Var_t*,const Args_t *>::template funBody<i>;

   template<fitnessFun f>
   using fFunBody = typename HeuPrivate::fFunArea_GA<void,const Var_t*,const Args_t *,Fitness_t*>::template funBody<f>;

   template<crossoverFun c>
   using cFunBody = typename HeuPrivate::cFunArea_GA<void,const Var_t*,const Var_t*,Var_t*,Var_t*,const Args_t *>::template funBody<c>;

   template<mutateFun m>
   using mFunBody= typename HeuPrivate::mFunArea_GA<void,Var_t *,const Args_t *>::template funBody<m>;

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
    using mutateFun = initializeFun;

    using ArgsType = void;

    template<initializeFun i>
    using iFunBody = typename HeuPrivate::iFunArea_GA<void,Var_t*>::template funBody<i>;

    template<fitnessFun f>
    using fFunBody = typename HeuPrivate::fFunArea_GA<void,const Var_t*,Fitness_t*>::template funBody<f>;

    template<crossoverFun c>
    using cFunBody = typename HeuPrivate::cFunArea_GA<void,const Var_t*,const Var_t*,Var_t*,Var_t*>::template funBody<c>;

    template<mutateFun m>
    using mFunBody= typename HeuPrivate::mFunArea_GA<void,Var_t *>::template funBody<m>;

    static const bool HasParameters=false;
};

#define Heu_MAKE_GAABSTRACT_TYPES \
using initializeFun = typename Base_t::initializeFun; \
using fitnessFun = typename Base_t::fitnessFun; \
using crossoverFun = typename Base_t::crossoverFun; \
using mutateFun = typename Base_t::mutateFun; \
using ArgsType = typename Base_t::ArgsType;


}   //  Heu

#endif  //  Heu_GAABSTRACT_HPP
