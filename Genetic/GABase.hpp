/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef GABASE_H
#define GABASE_H

#include "./GAOption.hpp"
#include <tuple>
#include <vector>
#include <list>
#include <cmath>
#include <random>
#include <algorithm>
#ifndef OptimT_NO_STATICASSERT
#include <type_traits>
#endif

#ifndef OptimT_NO_OUTPUT
#include <iostream>
#endif

#include "../OptimTemplates/Global"

namespace OptimT {

///Genetic algorithm base class
template<typename Var_t,typename Fitness_t,RecordOption Record,class ...Args>
class GABase {
public:
    ///Use std::tuple to store all extra parameters such as training samples
    using ArgsType = std::tuple<Args...>;
    ///Function to initialize Var
    using initializeFun = void(*)(Var_t*,const ArgsType*);
    ///Function to calculate fitness for Var
    using fitnessFun = void(*)(const Var_t*,const ArgsType*,Fitness_t*);
    ///Function to apply crossover for Var
    using crossoverFun = void(*)(const Var_t*,const Var_t*,Var_t*,Var_t*,const ArgsType*);
    ///Function to apply mutate for Var
    using mutateFun = initializeFun;

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


    ///Function to modify Args after each generation
    using otherOptFun = void(*)
        (ArgsType*,std::list<Gene>*,size_t generation,size_t failTimes,const GAOption*);
    ///list iterator to Gene
    using GeneIt_t = typename std::list<Gene>::iterator;

public:
    GABase() {

    };
    virtual ~GABase() {};

    ///initialize with options, initializeFun, fitnessFun, crossoverFun, mutateFun and Args
    virtual void initialize(
                            initializeFun _iFun,
                            fitnessFun _fFun,
                            crossoverFun _cFun,
                            mutateFun _mFun,
                            otherOptFun _ooF=nullptr,
                            const GAOption & options=GAOption(),
                            const ArgsType & args=ArgsType()) {
        _option=options;
        _population.resize(_option.populationSize);
        _args=args;
        _initializeFun=_iFun;
        _fitnessFun=_fFun;
        _crossoverFun=_cFun;
        _mutateFun=_mFun;

        if(_ooF==nullptr) {
            _otherOptFun=[](ArgsType*,std::list<Gene>*,size_t,size_t,const GAOption*){};
        } else {
            _otherOptFun=_ooF;
        }

        for(auto & i : _population) {
            _initializeFun(&i.self,&_args);
            i.setUncalculated();
        }
        customOptWhenInitialization();
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
#ifndef OptimT_NO_OUTPUT
                std::cout<<"Terminated by max generation limit"<<std::endl;
#endif
                break;
            }
            if(_option.maxFailTimes>0&&_failTimes>_option.maxFailTimes) {
#ifndef OptimT_NO_OUTPUT
                std::cout<<"Terminated by max failTime limit"<<std::endl;
#endif
                break;
            }
#ifndef OptimT_NO_OUTPUT
            std::cout<<"Generation "<<_generation
                    //<<" , elite fitness="<<_eliteIt->fitness()
                   <<std::endl;
#endif
            _otherOptFun(&_args,&_population,_generation,_failTimes,&_option);
            crossover();
            mutate();
        }
        _generation--;
    }

    ///best fitness
    //virtual Fitness_t bestFitness() const=0;

    ///the whole population
    const std::list<Gene> & population() const {
        return _population;
    }
    ///get option
    const GAOption & option() const {
        return _option;
    }
    ///generations used
    size_t generation() const {
        return _generation;
    }
    ///fail times
    size_t failTimes() const {
        return _failTimes;
    }
    ///other parameters
    const ArgsType & args() const {
        return _args;
    }

protected:
    std::list<Gene> _population;
    GAOption _option;

    size_t _generation;
    size_t _failTimes;

    ArgsType _args;

    fitnessFun _fitnessFun;
    initializeFun _initializeFun;
    crossoverFun _crossoverFun;
    mutateFun _mutateFun;
    otherOptFun _otherOptFun;
    virtual void customOptWhenInitialization() {}

    virtual void calculateAll() {
#ifdef OptimT_DO_PARALLELIZE
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
                _fitnessFun(&ptr->self,&_args,&ptr->_Fitness);
                ptr->_isCalculated=true;
            }
        }
#else
        for(Gene & i : _population) {
            if(i._isCalculated) {
                continue;
            }
            _fitnessFun(&i.self,&_args,&i._Fitness);
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

        while(crossoverQueue.empty()) {
            GeneIt_t a,b;
            a=crossoverQueue.back();
            crossoverQueue.pop_back();
            b=crossoverQueue.back();
            crossoverQueue.pop_back();
            _population.emplace_back();
            Gene * childA=&_population.back();
            _population.emplace_back();
            Gene * childB=&_population.back();
            _crossoverFun(&a->self,&b->self,&childA->self,&childB->self,&_args);
            childA->setUncalculated();
            childB->setUncalculated();
        }
    }


    ///mutate
    virtual void mutate()=0;


};



#define OPTIMT_MAKE_GABASE_TYPES \
using Gene = typename Base_t::Gene; \
using GeneIt_t = typename Base_t::GeneIt_t ; \
using initializeFun = typename Base_t::initializeFun; \
using fitnessFun = typename Base_t::fitnessFun; \
using crossoverFun = typename Base_t::crossoverFun; \
using mutateFun = typename Base_t::mutateFun; \
using otherOptFun = typename Base_t::otherOptFun; \
using ArgsType = typename Base_t::ArgsType; 



///partial specialization for GABase with record
template<typename Var_t,typename Fitness_t,class ...Args>
class GABase<Var_t,Fitness_t,RECORD_FITNESS,Args...>
    : public GABase<Var_t,Fitness_t,DONT_RECORD_FITNESS,Args...>
{
public:
    using Base_t = GABase<Var_t,Fitness_t,DONT_RECORD_FITNESS,Args...>;
    OPTIMT_MAKE_GABASE_TYPES

public:
    GABase() {

    };
    virtual ~GABase() {};
    ///initialize with options, initializeFun, fitnessFun, crossoverFun, mutateFun and Args
    virtual void initialize(
                            initializeFun _iFun,
                            fitnessFun _fFun,
                            crossoverFun _cFun,
                            mutateFun _mFun,
                            otherOptFun _ooF=nullptr,
                            const GAOption & options=GAOption(),
                            const ArgsType & args=ArgsType()) {
        this->_option=options;
        this->_population.resize(this->_option.populationSize);
        this->_args=args;
        this->_initializeFun=_iFun;
        this->_fitnessFun=_fFun;
        this->_crossoverFun=_cFun;
        this->_mutateFun=_mFun;

        if(_ooF==nullptr) {
           this-> _otherOptFun=[](ArgsType*,std::list<Gene>*,size_t,size_t,const GAOption*){};
        } else {
            this->_otherOptFun=_ooF;
        }

        for(auto & i : this->_population) {
            this->_initializeFun(&i.self,&this->_args);
            i.setUncalculated();
        }
        //_eliteIt=_population.begin();
        _record.clear();
        _record.reserve(this->_option.maxGenerations+1);

        customOptWhenInitialization();
    }


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
#ifndef OptimT_NO_OUTPUT
                std::cout<<"Terminated by max generation limit"<<std::endl;
#endif
                break;
            }
            if(this->_option.maxFailTimes>0&&this->_failTimes>this->_option.maxFailTimes) {
#ifndef OptimT_NO_OUTPUT
                std::cout<<"Terminated by max failTime limit"<<std::endl;
#endif
                break;
            }
#ifndef OptimT_NO_OUTPUT
            std::cout<<"Generation "<<this->_generation
                    //<<" , elite fitness="<<_eliteIt->fitness()
                   <<std::endl;
#endif
            this->_otherOptFun(&this->_args,
                               &this->_population,
                               this->_generation,
                               this->_failTimes,
                               &this->_option);
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

    virtual void customOptWhenInitialization() {}

};

}

#endif // GABASE_H
