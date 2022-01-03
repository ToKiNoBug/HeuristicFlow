/*
 Copyright © 2021  TokiNoBug
This file is part of AlgoTemplates.

    AlgoTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AlgoTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AlgoTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef GENETIC_H
#define GENETIC_H
#include <stdint.h>
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>

#ifdef OUTPUT
#include <iostream>
#endif
namespace AT {

};

namespace AT
{
///uniform random number in range [0,1]
double randD() {
    return (std::rand()%RAND_MAX)/double(RAND_MAX);
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

///Genetic algorithm
template<typename Var,bool isGreaterBetter,class ...Args>
class GA {
public:
    ///Use std::tuple to store all extra parameters such as training samples
    typedef std::tuple<Args...> ArgsType;
    ///Function to initialize Var
    typedef void(* initializeFun)(Var*,const ArgsType*);
    ///Function to calculate fitness for Var
    typedef double(* fitnessFun)(const Var*,const ArgsType*);
    ///Function to apply crossover for Var
    typedef void(* crossoverFun)(Var*,Var*,const ArgsType*);
    ///Function to apply mutate for Var
    typedef initializeFun mutateFun;

    ///Gene type for Var
    class Gene {
    public:
        Var self;
        bool isCalculated() const {
            return _isCalculated;
        }
        void setUncalculated() {
            _isCalculated=false;
        }
        double fitness() const {
            return _Fitness;
        }
    protected:
        friend class GA;
        bool _isCalculated;
        double _Fitness;
    };

    ///Function to modify Args after each generation
    typedef void(* otherOptFun)
        (ArgsType*,std::vector<Gene>*,size_t generation,size_t failTimes,const GAOption*);

public:
    GA() {
        _initializeFun=[](Var*,const ArgsType*){};
        _fitnessFun=[](const Var*,const ArgsType*){return 0.0;};
        _crossoverFun=[](Var*,Var*,const ArgsType*){};
        _mutateFun=[](Var*,const ArgsType*){};
        _otherOptFun=[](ArgsType*,std::vector<Gene>*,size_t,size_t,const GAOption*){};
    };
    virtual ~GA() {};
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
            _otherOptFun=[](ArgsType*,std::vector<Gene>*,size_t,size_t,const GAOption*){};
        } else {
            _otherOptFun=_ooF;
        }

        for(auto & i : _population) {
            _initializeFun(&i.self,&_args);
        }
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
#ifdef OUTPUT
                std::cout<<"Terminate by max generation limit"<<std::endl;
#endif
                break;
            }
            if(_option.maxFailTimes>0&&_failTimes>_option.maxFailTimes) {
#ifdef OUTPUT
                std::cout<<"Terminate by max failTime limit"<<std::endl;
#endif
                break;
            }
#ifdef OUTPUT
            std::cout<<"Generation "<<_generation<<std::endl;
#endif
            _otherOptFun(&_args,&_population,_generation,_failTimes,&_option);
            crossover();
            mutate();
        }
    }
    ///index of the best gene
    size_t eliteIdx() const {
        return _eliteIdx;
    }
    ///result
    const Var & result() const {
        return _population[_eliteIdx].self;
    }
    ///the whole population
    const std::vector<Gene> & population() const {
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
    size_t _eliteIdx;
    std::vector<Gene> _population;
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

    static inline bool isBetter(double A,double B) {
        if(isGreaterBetter) {
            return A>B;
        }
        return A<B;
    }

    virtual void select() {
        size_t bestIdx=0,worstIdx=0;
        double bestFitness=_population.front()._Fitness,worstFitness=bestFitness;
        for(size_t idx=1;idx<_population.size();idx++) {
            double curFitness=_population[idx]._Fitness;
            if(isBetter(curFitness,bestFitness)) {
                bestIdx=idx;
                bestFitness=curFitness;
            }
            if(isBetter(worstFitness,curFitness)) {
                worstIdx=idx;
                worstFitness=curFitness;
            }
        }
        if(_eliteIdx==bestIdx) {
            _failTimes++;
        } else {
            _failTimes=0;
            _eliteIdx=bestIdx;
        }

        _population[worstIdx]=_population[_eliteIdx];

    }
    virtual void crossover() {
        std::vector<Gene*> crossoverQueue;
        crossoverQueue.clear();
        crossoverQueue.reserve(_population.size());
        for(size_t idx=0;idx<_population.size();idx++) {
            if(idx==_eliteIdx) {
                continue;
            }

            if(randD()<=_option.crossoverProb) {
                crossoverQueue.push_back(idx+_population.data());
            }
        }

        std::random_shuffle(crossoverQueue.begin(),crossoverQueue.end());

        if(crossoverQueue.size()%2==1) {
            crossoverQueue.pop_back();
        }

        while(crossoverQueue.empty()) {
            Gene* a,*b;
            a=crossoverQueue.back();
            crossoverQueue.pop_back();
            b=crossoverQueue.back();
            crossoverQueue.pop_back();
            _crossoverFun(&a->self,&b->self,&_args);
            a->setUncalculated();
            b->setUncalculated();
        }
    }
    ///mutate
    virtual void mutate() {
        for(size_t idx=0;idx<_population.size();idx++) {
            if(idx==_eliteIdx) {
                continue;
            }
            if(randD()<=_option.mutateProb) {
                _mutateFun(&_population[idx].self,&_args);
                _population[idx].setUncalculated();
            }
        }
    }


};
}
#endif // GENETIC_H
