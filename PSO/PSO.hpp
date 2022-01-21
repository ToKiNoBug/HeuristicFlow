#ifndef PSO_HPP
#define PSO_HPP
#include "PSOOption.hpp"
#include <array>
#include <vector>
#include <tuple>
#ifdef EIGEN_CORE_H
#define OptimT_USER_USE_EIGEN
#endif

namespace OptimT {

template<size_t Dim,class...Args>
class PSO
{
public:
    using Var_t = std::array<double,Dim>;
    using Args_t = std::tuple<Args...>;
    using Fitness_t = double;
    using iFun_t = void(*)(Var_t * ,Var_t * ,const Args_t *);
    using fFun_t = void(*)(const Var_t * ,const Args_t *,Fitness_t *);

    class Particle
    {
    public:
        Var_t position;
        Var_t velocity;
        bool isCalculated;
        Fitness_t fitness;
    };
    using ooFun_t = void(*)
        (Args_t*,std::vector<Particle>*,size_t generation,size_t failTimes,PSOOption*);

public:
    PSO() {

    }

    virtual ~PSO() {};


    const PSOOption & option() {
        return _option;
    }

    virtual void initialize(iFun_t _i,
                            fFun_t _f,
                            ooFun_t _oo=default_ooFun,
                            const PSOOption & opt=PSOOption(),
                            const Args_t & args=Args_t()) {
        _iFun=_i;
        _fFun=_f;
        _ooFun=_oo;
        _option=opt;
        _args=args;
        _population.resize(_option.populationSize);

        for(Particle & i : _population) {
            _iFun(&i.position,&i.velocity,&_args);
        }
        _generation=0;
        _failTimes=0;
    }

    static void default_ooFun(Args_t*,
                              std::vector<Particle>*,
                              size_t generation,
                              size_t failTimes,
                              PSOOption*) {}

protected:
    PSOOption _option;
    Args_t _args;
    size_t _generation;
    size_t _failTimes;
    Var_t _posMin;
    Var_t _posMax;
    Var_t _velocityMax;

    iFun_t _iFun;
    fFun_t _fFun;
    ooFun_t _ooFun;



    std::vector<Particle> _population;



};

}

#endif // PSO_HPP
