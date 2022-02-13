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

#ifndef Heu_GABASE_H
#define Heu_GABASE_H
#include <iostream>
#include "./GAOption.hpp"
#include <tuple>
#include <vector>
#include <list>
#include <cmath>
#include <random>
#include <algorithm>
#include "GAAbstract.hpp"
#ifndef Heu_NO_STATICASSERT
#include <type_traits>
#endif

#ifdef Heu_DO_OUTPUT
#include <iostream>
#endif

#include "../HeuristicFlow/Global"

namespace Heu {


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
    RecordOption Record=DONT_RECORD_FITNESS,
    class Args_t=void>
class GABase : public GAAbstract<Var_t,Fitness_t,Args_t>
{
public:

    using Base_t = GAAbstract<Var_t,Fitness_t,Args_t>;

    Heu_MAKE_GAABSTRACT_TYPES

    ///Gene type for Var
    class Gene {
    public:
        Var_t self;
        bool isCalculated() const {
            return _isCalculated;
        }
        void setUncalculated() {
            _isCalculated=false;
        }
        Fitness_t fitness() const {
            return _Fitness;
        }

        bool _isCalculated;
        Fitness_t _Fitness;
    };
    ///list iterator to Gene
    using GeneIt_t = typename std::list<Gene>::iterator;

public:
    GABase() {
        _iFun=nullptr;
        _fFun=nullptr;
        _fFun=nullptr;
        _mFun=nullptr;
    };
    virtual ~GABase() {};

    ///initialize with options, initializeFun, fitnessFun, crossoverFun, mutateFun and Args

    inline void setiFun(initializeFun i) {
        _iFun=i;
    }

    inline initializeFun iFun() const {
        return _iFun;
    }

    inline void setfFun(fitnessFun f) {
        _fFun=f;
    }

    inline fitnessFun fFun() const {
        return _fFun;
    }
    
    inline void setcFun(crossoverFun c) {
        _cFun=c;
    }

    inline crossoverFun cFun() const {
        return _cFun;
    }

    inline void setmFun(mutateFun m) {
        _mFun=m;
    }

    inline mutateFun mFun() const {
        return _mFun;
    }

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

            if constexpr (Base_t::HasParameters)
                _iFun(&i.self,&this->_args);
            else
                _iFun(&i.self);

            i.setUncalculated();
        }
        customOptAfterInitialization();
    }
    
    ///start to solve
    virtual void run() {
        _generation=0;
        _failTimes=0;
        while(true) {
            _generation++;
            calculateAll();
            select();

            if(_generation>_option.maxGenerations) {
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
            customOptAfterEachGeneration();
            crossover();
            mutate();
        }
        _generation--;
    }

    ///best fitness
    //virtual Fitness_t bestFitness() const=0;

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

    fitnessFun _fFun;
    initializeFun _iFun;
    crossoverFun _cFun;
    mutateFun _mFun;
    virtual void customOptAfterInitialization() {}
    virtual void customOptAfterEachGeneration() {}

    virtual void calculateAll() {
#ifdef Heu_USE_THREADS
        static const uint32_t thN=OtGlobal::threadNum();
        std::vector<Gene*> tasks;
        tasks.resize(0);
        tasks.reserve(_population.size());
        for(Gene & i : _population) {
            if(i._isCalculated) {
                continue;
            }
            tasks.emplace_back(&i);
        }
#pragma omp parallel for
        for(uint32_t begIdx=0;begIdx<thN;begIdx++) {
            for(uint32_t i=begIdx;i<tasks.size();i+=thN) {
                Gene * ptr=tasks[i];
                if constexpr (Base_t::HasParameters)
                    _fFun(&ptr->self,&this->_args,&ptr->_Fitness);
                else {
                    _fFun(&ptr->self,&ptr->_Fitness);
                }
                ptr->_isCalculated=true;
            }
        }
#else
        for(Gene & i : _population) {
            if(i._isCalculated) {
                continue;
            }
            if constexpr (Base_t::HasParameters)
                _fFun(&i.self,&this->_args,&i._Fitness);
            else
                _fFun(&i.self,&i._Fitness);

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
                    global_mt19937);

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
            if constexpr (Base_t::HasParameters)
                _cFun(&a->self,&b->self,&childA->self,&childB->self,&this->_args);
            else
                _cFun(&a->self,&b->self,&childA->self,&childB->self);

            childA->setUncalculated();
            childB->setUncalculated();
        }
    }


    ///mutate
    virtual void mutate()=0;


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
template<typename Var_t,typename Fitness_t,class Args_t>
class GABase<Var_t,Fitness_t,RECORD_FITNESS,Args_t>
    : public GABase<Var_t,Fitness_t,DONT_RECORD_FITNESS,Args_t>
{
public:
    using Base_t = GABase<Var_t,Fitness_t,DONT_RECORD_FITNESS,Args_t>;
    Heu_MAKE_GABASE_TYPES

public:
    GABase() {};

    virtual ~GABase() {};

    ///start to solve
    virtual void run() {
        this->_generation=0;
        this->_failTimes=0;


        _record.clear();
        _record.reserve(this->_option.maxGenerations+1);

        while(true) {
            this->_generation++;
            this->calculateAll();
            this->select();
            _record.emplace_back(bestFitness());
            if(this->_generation>this->_option.maxGenerations) {
#ifdef Heu_DO_OUTPUT
                std::cout<<"Terminated by max generation limit"<<std::endl;
#endif
                break;
            }
            if(this->_option.maxFailTimes>0&&this->_failTimes>this->_option.maxFailTimes) {
#ifdef Heu_DO_OUTPUT
                std::cout<<"Terminated by max failTime limit"<<std::endl;
#endif
                break;
            }
#ifdef Heu_DO_OUTPUT
            std::cout<<"Generation "<<this->_generation
                    //<<" , elite fitness="<<_eliteIt->fitness()
                   <<std::endl;
#endif      
            this->customOptAfterEachGeneration();
            this->crossover();
            this->mutate();
        }
        this->_generation--;
    }

    ///best fitness
    virtual Fitness_t bestFitness() const=0;

    ///result
    const std::vector<Fitness_t> record() const {
        return _record;
    }

protected:
    std::vector<Fitness_t> _record;

};

}

#endif // GABASE_H
