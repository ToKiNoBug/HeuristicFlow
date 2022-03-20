// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_GABASE_H
#define EIGEN_HEU_GABASE_H
#include "./GAOption.hpp"
#include <tuple>
#include <vector>
#include <list>
#include <cmath>
#include <random>
#include <algorithm>
#include "GAAbstract.hpp"
#include <type_traits>

#ifdef EIGEN_HEU_DO_OUTPUT
#include <iostream>
#endif

#include "../../Global"

namespace Eigen
{

namespace internal
{
    /**
   *  @brief Genetic algorithm base class. 
   *  It's an abstrcat base class for all genetic algorithm solvers.
   *
   *  @author TokiNoBug
   * 
   *  @tparam Var_t  Type of decisition variable
   *  @tparam Fitness_t  Type of fitness value(objective value)
   *  @tparam Record  Whether the solver records fitness changelog
   *  @tparam ...Args  Type of other parameters.
   * 
  */
template<typename Var_t,typename Fitness_t,
    RecordOption Record,
    class Args_t,
         typename GAAbstract<Var_t,Fitness_t,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,Fitness_t,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,Fitness_t,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,Fitness_t,Args_t>::mutateFun _mFun_>
class GABase :
    public GAAbstract<Var_t,Fitness_t,Args_t> ,
    public GAAbstract<Var_t,Fitness_t,Args_t>::template iFunBody<_iFun_>,
    public GAAbstract<Var_t,Fitness_t,Args_t>::template fFunBody<_fFun_>,
    public GAAbstract<Var_t,Fitness_t,Args_t>::template cFunBody<_cFun_>,
    public GAAbstract<Var_t,Fitness_t,Args_t>::template mFunBody<_mFun_>
{
private:
    using Base_t = GAAbstract<Var_t,Fitness_t,Args_t>;
public:
    Heu_MAKE_GAABSTRACT_TYPES

    ///Gene type for Var
    class Gene {
    public:
    using fastFitness_t = typename
        std::conditional<(sizeof(Fitness_t)>sizeof(double)),const Fitness_t &,Fitness_t>::type;
        Var_t self;
        bool isCalculated() const {
            return _isCalculated;
        }
        void setUncalculated() {
            _isCalculated=false;
        }
        fastFitness_t fitness() const {
            return _Fitness;
        }

        bool _isCalculated;
        Fitness_t _Fitness;
    };
    ///list iterator to Gene
    using GeneIt_t = typename std::list<Gene>::iterator;

public:
    GABase() {
    };
    virtual ~GABase() {};

    ///initialize with options, initializeFun, fitnessFun, crossoverFun, mutateFun and Args

    inline void setOption(const GAOption & o) {
        _option=o;
    }

    ///get option
    inline const GAOption & option() const {
        return _option;
    }
    
    void initializePop() {
        _population.resize(_option.populationSize);
        for(auto & i : _population) {
            GAExecutor<Base_t::HasParameters>::doInitialization(this,&i.self);

            i.setUncalculated();
        }
        customOptAfterInitialization();
    }
    
    ///start to solve
    template<class this_t=GABase>
    void run() {
        _generation=0;
        _failTimes=0;
        static_cast<this_t*>(this)->__impl_clearRecord();
        while(true) {
            _generation++;
            calculateAll();

            select();

            static_cast<this_t*>(this)->__impl_recordFitness();

            if(_generation>_option.maxGenerations) {
#ifdef EIGEN_HEU_DO_OUTPUT
                std::cout<<"Terminated by max generation limitation"<<std::endl;
#endif
                break;
            }
            if(_option.maxFailTimes>0&&_failTimes>_option.maxFailTimes) {
#ifdef EIGEN_HEU_DO_OUTPUT
                std::cout<<"Terminated by max failTime limitation"<<std::endl;
#endif
                break;
            }
#ifdef EIGEN_HEU_DO_OUTPUT
            std::cout<<"Generation "<<_generation
                    //<<" , elite fitness="<<_eliteIt->fitness()
                   <<std::endl;
#endif
            customOptAfterEachGeneration();
            crossover();
            mutate();
        }
        _generation--;
    }

    ///the whole population
    inline const std::list<Gene> & population() const {
        return _population;
    }
    ///generations used
    inline size_t generation() const {
        return _generation;
    }
    ///fail times
    inline size_t failTimes() const {
        return _failTimes;
    }

protected:
    std::list<Gene> _population;
    GAOption _option;

    size_t _generation;
    size_t _failTimes;

    virtual void customOptAfterInitialization() {}
    virtual void customOptAfterEachGeneration() {}

    inline void __impl_clearRecord() {}
    inline void __impl_recordFitness() {}

    virtual void calculateAll() {
#ifdef EIGEN_HAS_OPENMP
        std::vector<Gene*> tasks;
        tasks.resize(0);
        tasks.reserve(_population.size());
        for(Gene & i : _population) {
            if(i._isCalculated) {
                continue;
            }
            tasks.emplace_back(&i);
        }
        static const int32_t thN=Eigen::nbThreads();
        #pragma omp parallel for schedule(dynamic,tasks.size()/thN)
        for(int i=0;i<tasks.size();i++) {
                Gene * ptr=tasks[i];

                GAExecutor<Base_t::HasParameters>::doFitness(this,&ptr->self,&ptr->_Fitness);

                ptr->_isCalculated=true;
        }
#else
        for(Gene & i : _population) {
            if(i._isCalculated) {
                continue;
            }

            GAExecutor<Base_t::HasParameters>::doFitness(this,&i.self,&i._Fitness);

            i._isCalculated=true;
        }
#endif
    }

    virtual void select()=0;

    virtual void crossover() {
        std::vector<GeneIt_t> crossoverQueue;
        crossoverQueue.clear();
        crossoverQueue.reserve(_population.size());

        for(GeneIt_t it=_population.begin();it!=_population.end();it++) {
            if(randD()<=_option.crossoverProb) {
                crossoverQueue.emplace_back(it);
            }
        }


        std::shuffle(
                    crossoverQueue.begin(),
                    crossoverQueue.end(),
                    global_mt19937());

        if(crossoverQueue.size()%2==1) {
            crossoverQueue.pop_back();
        }

        while(!crossoverQueue.empty()) {
            GeneIt_t a,b;
            a=crossoverQueue.back();
            crossoverQueue.pop_back();
            b=crossoverQueue.back();
            crossoverQueue.pop_back();
            _population.emplace_back();
            Gene * childA=&_population.back();
            _population.emplace_back();
            Gene * childB=&_population.back();

            GAExecutor<Base_t::HasParameters>::
                    doCrossover(this,&a->self,&b->self,&childA->self,&childB->self);

            childA->setUncalculated();
            childB->setUncalculated();
        }
    }


    ///mutate
    virtual void mutate() {
        std::vector<GeneIt_t> mutateList;
        mutateList.reserve(size_t(this->_population.size()*this->_option.mutateProb*2));
        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            if(randD()<=this->_option.mutateProb) {
                mutateList.emplace_back(it);
            }
        }
        for(auto src : mutateList) {
            this->_population.emplace_back();
            GAExecutor<Base_t::HasParameters>::doMutation(this,&src->self,&this->_population.back().self);
            this->_population.back().setUncalculated();
        }
    }

protected:
    template<bool HasParameters,class unused=void>
    struct GAExecutor
    {
        inline static void doInitialization(GABase * s,Var_t * v) {
            s->runiFun(v,&s->_args);
        }

        inline static void doFitness(GABase * s,const Var_t * v,Fitness_t * f) {
            s->runfFun(v,&s->_args,f);
        }

        inline static void doCrossover(GABase * s,const Var_t * p1,const Var_t *p2,
                                Var_t *c1,Var_t * c2) {
            s->runcFun(p1,p2,c1,c2,&s->_args);
        }

        inline static void doMutation(GABase * s,const Var_t * src,Var_t * dst) {
            s->runmFun(src,dst,&s->_args);
        }
        static_assert (HasParameters==GABase::HasParameters,
                "struct GAExecutor actived with wrong template parameter.");
    };

    template<class unused>
    struct GAExecutor<false,unused>
    {
        inline static void doInitialization(GABase * s,Var_t * v) {
            s->runiFun(v);
        }

        inline static void doFitness(GABase * s,const Var_t * v,Fitness_t * f) {
            s->runfFun(v,f);
        }

        inline static void doCrossover(GABase * s,const Var_t * p1,const Var_t *p2,
                                Var_t *c1,Var_t * c2) {
            s->runcFun(p1,p2,c1,c2);
        }

        inline static void doMutation(GABase * s,const Var_t * src,Var_t * dst) {
            s->runmFun(src,dst);
        }

        static_assert(GABase::HasParameters==false,
                "struct GAExecutor actived with wrong template parameter.");
    };



};

#define Heu_MAKE_GABASE_TYPES \
using Gene = typename Base_t::Gene; \
using GeneIt_t = typename Base_t::GeneIt_t; \
Heu_MAKE_GAABSTRACT_TYPES

/**
   *  @brief partial specialization for GABase with record.
   *  It's an abstrcat base class for all genetic algorithm solvers that records fitness.
   *
   *  @tparam Var_t  Type of decisition variable
   *  @tparam Fitness_t  Type of fitness value(objective value)
   *  @tparam RecordOption  Whether the solver records fitness changelog
   *  @tparam Args_t  Type of other parameters.
  */
template<typename Var_t,typename Fitness_t,class Args_t,
         typename GAAbstract<Var_t,Fitness_t,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,Fitness_t,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,Fitness_t,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,Fitness_t,Args_t>::mutateFun _mFun_>
class GABase<Var_t,Fitness_t,RECORD_FITNESS,Args_t,_iFun_,_fFun_,_cFun_,_mFun_>
    : public GABase<Var_t,Fitness_t,DONT_RECORD_FITNESS,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>
{
private:
    using Base_t = GABase<Var_t,Fitness_t,DONT_RECORD_FITNESS,Args_t,
                                            _iFun_,_fFun_,_cFun_,_mFun_>;
    friend Base_t;
public:
    Heu_MAKE_GABASE_TYPES
    GABase() {};

    virtual ~GABase() {};

    ///start to solve
    template<class this_t=GABase>
    inline void run() {
        Base_t::template run<this_t>();
    }

    ///best fitness
    virtual Fitness_t bestFitness() const=0;

    ///result
    const std::vector<Fitness_t> record() const {
        return _record;
    }


protected:
    std::vector<Fitness_t> _record;

    inline void __impl_clearRecord() {
        _record.clear();
        _record.reserve(this->_option.maxGenerations+1);
    }

    inline void __impl_recordFitness() {
        _record.emplace_back(bestFitness());
    }

};
}   //  internal
}   //  namespace Eigen

#endif // EIGEN_HEU_GABASE_H
