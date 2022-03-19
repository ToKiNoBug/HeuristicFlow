// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_PSOABSTRCAT_HPP
#define Heu_PSOABSTRCAT_HPP

#include "../../Global"
#include "PSOOption.hpp"
#include "PSOParameterPack.hpp"

#ifdef Heu_DO_OUTPUT
#include <iostream>
#endif

namespace Heu
{

namespace internal
{
    
/**
 * @brief Base-base class for PSO solvers. Some fundamental typedefs and functions here.
 * 
 * @tparam Var_t Type of determination vector.
 * @tparam Record Record trainning curve or not.
 * @tparam Arg_t Any other parameters
 */
template<class Var_t,
    class Fitness_t,
    RecordOption Record,
    class Arg_t,
         typename PSOParameterPack<Var_t,Fitness_t,Arg_t>::iFun_t _iFun_,
         typename PSOParameterPack<Var_t,Fitness_t,Arg_t>::fFun_t _fFun_>
class PSOAbstract :
        public PSOParameterPack<Var_t,Fitness_t,Arg_t>,
        public PSOParameterPack<Var_t,Fitness_t,Arg_t>::template iFunBody<_iFun_>,
        public PSOParameterPack<Var_t,Fitness_t,Arg_t>::template fFunBody<_fFun_>
{
    using Base_t = PSOParameterPack<Var_t,Fitness_t,Arg_t>;
public:
    Heu_MAKE_PSOPARAMETERPACK_TYPES
    
    class Point
    {
    public:
        Var_t position;
        Fitness_t fitness;
    };

    class Particle : public Point
    {
    public:
        Var_t velocity;
        //bool isCalculated;
        Point pBest;

    };

public:
    PSOAbstract() {};
    virtual ~PSOAbstract() {};

    inline void setOption(const PSOOption & opt) {
        _option=opt;
    }

    inline const PSOOption & option() const {
        return _option;
    }

    inline size_t generation() const {
        return _generation;
    }

    inline size_t failTimes() const {
        return _failTimes;
    }

    inline const Var_t & posMin() const {
        return _posMin;
    }

    inline const Var_t & posMax() const {
        return _posMax;
    }

    inline const Var_t & velocityMax() const {
        return _velocityMax;
    }

    inline const std::vector<Particle> & population() const {
        return _population;
    }

    inline const Point & globalBest() const {
        return gBest;
    }

    inline void setPVRange(const Var_t & pMin,
            const Var_t & pMax,
            const Var_t & vMax) {
        _posMin=pMin;
        _posMax=pMax;
        _velocityMax=vMax;
    }

    void initializePop() {
        _population.resize(_option.populationSize);

        for(Particle & i : _population) {
            PSOExecutor<Base_t::HasParameters>::doInitialize(this,&i.position,
                                                             &i.velocity,&_posMin,
                                                             &_posMax,&_velocityMax);

            PSOExecutor<Base_t::HasParameters>::doFitness(this,&i.position,&i.fitness);

            i.pBest=i;
        }

        gBest=_population.front();
        _generation=0;
        _failTimes=0;
    }

    template<class this_t=PSOAbstract>
    void run() {
        _generation=0;
        _failTimes=0;

        static_cast<this_t*>(this)->__impl_clearRecord();

        while(true) {
            _generation++;
            calculateAll();
            updatePGBest();

            static_cast<this_t*>(this)->__impl_recordFitness();
            if(_generation>_option.maxGeneration) {
#ifdef Heu_DO_OUTPUT
                    std::cout<<"Terminated by max generation limit"<<std::endl;
#endif
                    break;
            }

            if(_option.maxFailTimes>0&&_failTimes>_option.maxFailTimes) {
#ifdef Heu_DO_OUTPUT
                    std::cout<<"Terminated by max failTime limit"<<std::endl;
#endif
                    break;
            }
#ifdef Heu_DO_OUTPUT
            std::cout<<"Generation "<<_generation
                        //<<" , elite fitness="<<_eliteIt->fitness()
                       <<std::endl;
#endif
            updatePopulation();

            customOptAfterEachGeneration();
        }
        _generation--;

    }

protected:
    PSOOption _option;
    size_t _generation;
    size_t _failTimes;

    Var_t _posMin;
    Var_t _posMax;
    Var_t _velocityMax;

    std::vector<Particle> _population;

    Point gBest;

    inline void __impl_clearRecord() {}

    inline void __impl_recordFitness() {}

    virtual void calculateAll() {
#ifdef Heu_USE_THREADS
        static const int32_t thN=Eigen::nbThreads();
#pragma omp parallel for schedule(dynamic,_population.size()/thN)
        for(int i=0;i<_population.size();i++) {
                Particle * ptr=&_population[i];
                PSOExecutor<Base_t::HasParameters>::
                        doFitness(this,&ptr->position,&ptr->fitness);
        }
#else
        for(Particle & i : _population) {
            PSOExecutor<Base_t::HasParameters>::doFitness(this,&i.position,&i.fitness);
        }
#endif
    }

    virtual void updatePGBest()=0;

    virtual void updatePopulation()=0;

    virtual void customOptAfterEachGeneration() {};


    template<bool _HasParameters,class unused=void>
    struct PSOExecutor
    {
        inline static void doInitialize(PSOAbstract * s,Var_t * pos,Var_t * velocity,
                                        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax) {
            s->runiFun(pos,velocity,pMin,pMax,vMax,&s->_arg);
        }

        inline static void doFitness(PSOAbstract * s,
                                     const Var_t * pos,Fitness_t * f) {
            s->runfFun(pos,&s->_arg,f);
        }

        static_assert (Base_t::HasParameters==_HasParameters,
            "A wrong specialization of PSOExecuter is called");
    };

    template<class unused>
    struct PSOExecutor<false,unused>
    {
        inline static void doInitialize(PSOAbstract * s,Var_t * pos,Var_t * velocity,
                                        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax) {
            s->runiFun(pos,velocity,pMin,pMax,vMax);
        }

        inline static void doFitness(PSOAbstract * s,
                                     const Var_t * pos,Fitness_t * f) {
            s->runfFun(pos,f);
        }

        static_assert (Base_t::HasParameters==false,
            "A wrong specialization of PSOExecuter is called");
    };
};

#define Heu_MAKE_PSOABSTRACT_TYPES \
Heu_MAKE_PSOPARAMETERPACK_TYPES \
using Point_t = typename Base_t::Point; \
using Particle_t = typename Base_t::Particle;





///partial specialization for PSO with recording
template<class Var_t,class Fitness_t,class Arg_t,
         typename PSOParameterPack<Var_t,Fitness_t,Arg_t>::iFun_t _iFun_,
         typename PSOParameterPack<Var_t,Fitness_t,Arg_t>::fFun_t _fFun_>
class PSOAbstract<Var_t,Fitness_t,RECORD_FITNESS,Arg_t,_iFun_,_fFun_>
        : public PSOAbstract<Var_t,Fitness_t,DONT_RECORD_FITNESS,Arg_t,_iFun_,_fFun_>
{
    using Base_t = PSOAbstract<Var_t,Fitness_t,DONT_RECORD_FITNESS,Arg_t,_iFun_,_fFun_>;
    friend Base_t;
public:
    Heu_MAKE_PSOABSTRACT_TYPES

    const std::vector<Fitness_t> & record() const {
        return _record;
    }

    virtual Fitness_t bestFitness() const=0;

    template<class this_t=PSOAbstract>
    void run() {
        Base_t::template run<this_t>();
    }

protected:
    std::vector<Fitness_t> _record;

    inline void __impl_clearRecord() {
        _record.clear();
        _record.reserve(this->_option.maxGeneration+1);
    }

    inline void __impl_recordFitness() {
        _record.emplace_back(bestFitness());
    }

};

}   //  internal

}   //  Heu


#endif // PSOABSTRCAT_HPP
