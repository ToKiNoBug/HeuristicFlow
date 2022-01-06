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

#include <stdint.h>
#include <cmath>
#include <random>
#include <list>
#include <algorithm>

#ifndef OptimT_NO_OUTPUT
#include <iostream>
#endif

namespace OptimT {

std::random_device global_random_device;
std::mt19937 global_mt19937(global_random_device());

const bool RECORD_FITNESS=true;
const bool DONT_RECORD_FITNESS=false;

///uniform random number in range [0,1)
double randD() {
    static std::uniform_real_distribution<double> rnd(0,1);
    static std::mt19937 mt(global_random_device());
    return rnd.operator()(mt);
    //return (std::rand()%RAND_MAX)/double(RAND_MAX);
}
///uniform random number in range [min,max]
double randD(const double min,const double max) {
    return (max-min)*randD()+min;
}

///options about GA algorithm
struct GAOption
{
public:
    GAOption() {
        populationSize=100;
        maxFailTimes=50;
        maxGenerations=300;
        crossoverProb=0.8;
        mutateProb=0.05;
    }
    size_t populationSize;
    size_t maxFailTimes;
    size_t maxGenerations;
    double crossoverProb;
    double mutateProb;
};

///Genetic algorithm base class
template<typename Var_t,typename Fitness_t,bool Record,class ...Args>
class GABase {
public:
    ///Use std::tuple to store all extra parameters such as training samples
    using ArgsType = std::tuple<Args...>;
    ///Function to initialize Var
    using initializeFun = void(*)(Var_t*,const ArgsType*);
    ///Function to calculate fitness for Var
    using fitnessFun = Fitness_t(*)(const Var_t*,const ArgsType*);
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

        _initializeFun=[](Var_t*,const ArgsType*){};
        _fitnessFun=[](const Var_t*,const ArgsType*){return 0.0;};
        _crossoverFun=[](const Var_t*,const Var_t*,Var_t*,Var_t*,const ArgsType*){};
        _mutateFun=[](Var_t*,const ArgsType*){};
        _otherOptFun=[](ArgsType*,std::list<Gene>*,size_t,size_t,const GAOption*){};

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
        _eliteIt=_population.begin();
    }


    ///start to solve
    virtual void run() {
        _generation=0;
        _failTimes=0;
        while(true) {
            _generation++;
            calculateAll();
            select();
            //_recording.emplace_back(_eliteIt->_Fitness);
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
            std::cout<<"Generation "<<_generation<<" , elite fitness="<<_eliteIt->fitness()<<std::endl;
#endif
            _otherOptFun(&_args,&_population,_generation,_failTimes,&_option);
            crossover();
            mutate();
        }
    }

    ///index of the best gene
    GeneIt_t eliteIt() const {
        return _eliteIt;
    }
    ///result
    const Var_t & result() const {
        return _eliteIt->self;
    }
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
    GeneIt_t _eliteIt;
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

    virtual void calculateAll() {
        for(Gene & i : _population) {
            if(i._isCalculated) {
                continue;
            }
            i._Fitness=_fitnessFun(&i.self,&_args);
            i._isCalculated=true;
        }
    }

    virtual void select()=0;

    virtual void crossover() {
        std::vector<GeneIt_t> crossoverQueue;
        crossoverQueue.clear();
        crossoverQueue.reserve(_population.size());

        for(GeneIt_t it=_population.begin();it!=_population.end();it++) {
            if(it==_eliteIt) {
                continue;
            }

            if(randD()<=_option.crossoverProb) {
                crossoverQueue.emplace_back(it);
            }
        }
        std::shuffle(crossoverQueue.begin(),crossoverQueue.end(),global_mt19937);

        if(crossoverQueue.size()%2==1) {
            crossoverQueue.pop_back();
        }

        while(crossoverQueue.empty()) {
            GeneIt_t a,b;
            a=crossoverQueue.back();
            crossoverQueue.pop_back();
            b=crossoverQueue.back();
            crossoverQueue.pop_back();
            Gene * childA=&_population.emplace_back();
            Gene * childB=&_population.emplace_back();
            _crossoverFun(&a->self,&b->self,&childA->self,&childB->self,&_args);
            childA->setUncalculated();
            childB->setUncalculated();
        }
    }


    ///mutate
    virtual void mutate() {
        for(auto it=_population.begin();it!=_population.end();++it) {
            if(it==_eliteIt) {
                continue;
            }
            if(randD()<=_option.mutateProb) {
                _mutateFun(&it->self,&_args);
                it->setUncalculated();
            }
        }
    }


};


///partial specialization for GABase with record
template<typename Var_t,typename Fitness_t,class ...Args>
class GABase<Var_t,Fitness_t,true,Args...>
{
public:
    ///Use std::tuple to store all extra parameters such as training samples
    using ArgsType = std::tuple<Args...>;
    ///Function to initialize Var
    using initializeFun = void(*)(Var_t*,const ArgsType*);
    ///Function to calculate fitness for Var
    using fitnessFun = Fitness_t(*)(const Var_t*,const ArgsType*);
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

        _initializeFun=[](Var_t*,const ArgsType*){};
        _fitnessFun=[](const Var_t*,const ArgsType*){return 0.0;};
        _crossoverFun=[](const Var_t*,const Var_t*,Var_t*,Var_t*,const ArgsType*){};
        _mutateFun=[](Var_t*,const ArgsType*){};
        _otherOptFun=[](ArgsType*,std::list<Gene>*,size_t,size_t,const GAOption*){};

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
        _eliteIt=_population.begin();
        _record.clear();
        _record.reserve(_option.maxGenerations+1);
    }


    ///start to solve
    virtual void run() {
        _generation=0;
        _failTimes=0;
        while(true) {
            _generation++;
            calculateAll();
            select();
            _record.emplace_back(_eliteIt->_Fitness);
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
            std::cout<<"Generation "<<_generation<<" , elite fitness="<<_eliteIt->fitness()<<std::endl;
#endif
            _otherOptFun(&_args,&_population,_generation,_failTimes,&_option);
            crossover();
            mutate();
        }
    }

    ///index of the best gene
    GeneIt_t eliteIt() const {
        return _eliteIt;
    }
    ///result
    const Var_t & result() const {
        return _eliteIt->self;
    }
    const std::vector<Fitness_t> record() const {
        return _record;
    }

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
    GeneIt_t _eliteIt;
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

    std::vector<Fitness_t> _record;

    virtual void calculateAll() {
        for(Gene & i : _population) {
            if(i._isCalculated) {
                continue;
            }
            i._Fitness=_fitnessFun(&i.self,&_args);
            i._isCalculated=true;
        }
    }

    virtual void select()=0;

    virtual void crossover() {
        std::vector<GeneIt_t> crossoverQueue;
        crossoverQueue.clear();
        crossoverQueue.reserve(_population.size());

        for(GeneIt_t it=_population.begin();it!=_population.end();it++) {
            if(it==_eliteIt) {
                continue;
            }

            if(randD()<=_option.crossoverProb) {
                crossoverQueue.emplace_back(it);
            }
        }

        std::shuffle(crossoverQueue.begin(),crossoverQueue.end(),global_mt19937);

        if(crossoverQueue.size()%2==1) {
            crossoverQueue.pop_back();
        }

        while(crossoverQueue.empty()) {
            GeneIt_t a,b;
            a=crossoverQueue.back();
            crossoverQueue.pop_back();
            b=crossoverQueue.back();
            crossoverQueue.pop_back();
            Gene * childA=&_population.emplace_back();
            Gene * childB=&_population.emplace_back();
            _crossoverFun(&a->self,&b->self,&childA->self,&childB->self,&_args);
            childA->setUncalculated();
            childB->setUncalculated();
        }
    }


    ///mutate
    virtual void mutate() {
        for(auto it=_population.begin();it!=_population.end();++it) {
            if(it==_eliteIt) {
                continue;
            }
            if(randD()<=_option.mutateProb) {
                _mutateFun(&it->self,&_args);
                it->setUncalculated();
            }
        }
    }


};


/*



template<typename Var_t,typename Fitness_t,class ...Args>
class GABase<Var_t,Fitness_t,true,Args...>
        : protected GABase<Var_t,Fitness_t,false,Args...>
{

private:
    using GABase<Var_t,Fitness_t,false,Args...> Base_t;

    ///Use std::tuple to store all extra parameters such as training samples
    using std::tuple<Args...> ArgsType;
    ///Function to initialize Var
    using void(* initializeFun)(Var_t*,const ArgsType*);
    ///Function to calculate fitness for Var
    using Fitness_t(* fitnessFun)(const Var_t*,const ArgsType*);
    ///Function to apply crossover for Var
    using void(* crossoverFun)(const Var_t*,const Var_t*,Var_t*,Var_t*,const ArgsType*);
    ///Function to apply mutate for Var
    using initializeFun mutateFun;
    ///Function to modify Args after each generation
    using Base_t::Gene Gene_t;
    using void(* otherOptFun)
        (ArgsType*,std::list<Gene_t>*,size_t generation,size_t failTimes,const GAOption*);
public:
    GABase() {};
    virtual ~GABase() {};

    const std::vector<Fitness_t> & recording() const {
        return _recording;
    }
    ///rewrite initialize function for doing record.
    virtual void initialize(
                            initializeFun _iFun,
                            fitnessFun _fFun,
                            crossoverFun _cFun,
                            mutateFun _mFun,
                            otherOptFun  _ooF=nullptr,
                            const GAOption & options=GAOption(),
                            const ArgsType & args=ArgsType()) {
        Base_t::_option=options;
        Base_t::_population.resize(Base_t::_option.populationSize);
        Base_t::_args=args;
        Base_t::_initializeFun=_iFun;
        Base_t::_fitnessFun=_fFun;
        Base_t::_crossoverFun=_cFun;
        Base_t::_mutateFun=_mFun;

        if(_ooF==nullptr) {
            Base_t::_otherOptFun=[](ArgsType*,typeof(Base_t::_population)*,size_t,size_t,const GAOption*){};
        } else {
            Base_t::_otherOptFun=_ooF;
        }

        for(auto & i : Base_t::_population) {
            _initializeFun(&i.self,&Base_t::_args);
            i.setUncalculated();
        }
        Base_t::_eliteIt=Base_t::_population.begin();
        _recording.clear();
        _recording.reserve(Base_t::_option.maxGenerations+1);
    }
    ///rewrite run function
    virtual void run() {
        Base_t::_generation=0;
        Base_t::_failTimes=0;
        while(true) {
            Base_t::_generation++;
            Base_t::calculateAll();
            Base_t::select();
            _recording.emplace_back(Base_t::_eliteIt->_Fitness);
            if(Base_t::_generation>Base_t::_option.maxGenerations) {
#ifndef OptimT_NO_OUTPUT
                std::cout<<"Terminated by max generation limit"<<std::endl;
#endif
                break;
            }
            if(Base_t::_option.maxFailTimes>0&&Base_t::_failTimes>Base_t::_option.maxFailTimes) {
#ifndef OptimT_NO_OUTPUT
                std::cout<<"Terminated by max failTime limit"<<std::endl;
#endif
                break;
            }
#ifndef OptimT_NO_OUTPUT
            std::cout<<"Generation "<<Base_t::_generation<<" , elite fitness="<<Base_t::_eliteIt->fitness()<<std::endl;
#endif
            _otherOptFun(&Base_t::_args,&Base_t::_population,Base_t::_generation,Base_t::_failTimes,&Base_t::_option);
            Base_t::crossover();
            Base_t::mutate();
        }
    }

protected:
    std::vector<Fitness_t> _recording;
};
*/
}

#endif // GABASE_H
