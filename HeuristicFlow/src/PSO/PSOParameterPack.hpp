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

    /*
    inline void doInitialize(iFun_t iFun,Var_t * pos,Var_t * velocity,
        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax) {
            iFun(pos,velocity,pMin,pMax,vMax,&_arg);
    }

    inline void doFitness(fFun_t fFun,const Var_t * v,Fitness_t * f) {
        fFun(v,&_arg,f);
    }
    */
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

    /*
    inline static void doInitialize(iFun_t iFun,Var_t * pos,Var_t * velocity,
        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax) {
            iFun(pos,velocity,pMin,pMax,vMax);
    }

    inline static void doFitness(fFun_t fFun,const Var_t * v,Fitness_t * f) {
        fFun(v,f);
    }
    */

    static const bool HasParameters=false;
};

#define Heu_MAKE_PSOPARAMETERPACK_TYPES \
using Args_t = typename Base_t::Args_t; \
using iFun_t = typename Base_t::iFun_t; \
using fFun_t = typename Base_t::fFun_t;

} // namespace Heu


#endif  //  Heu_PSOPARAMETERPACK_HPP
