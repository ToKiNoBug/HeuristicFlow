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

#ifndef GENETIC_H
#define GENETIC_H
#include <stdint.h>
#include <vector>
#include <list>
#include <tuple>
#include <cmath>
#include <random>
#include <algorithm>

#ifndef OptimT_NO_OUTPUT
#include <iostream>
#endif
namespace OptimT {

};

namespace OptimT
{
///uniform random number in range [0,1)
double randD() {
    static std::uniform_real_distribution<double> rnd(0,1);
    static std::mt19937 mt(std::rand());
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
    typedef void(* crossoverFun)(const Var*,const Var*,Var*,Var*,const ArgsType*);
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
        (ArgsType*,std::list<Gene>*,size_t generation,size_t failTimes,const GAOption*);

public:
    GA() {
        _initializeFun=[](Var*,const ArgsType*){};
        _fitnessFun=[](const Var*,const ArgsType*){return 0.0;};
        _crossoverFun=[](const Var*,const Var*,Var*,Var*,const ArgsType*){};
        _mutateFun=[](Var*,const ArgsType*){};
        _otherOptFun=[](ArgsType*,std::list<Gene>*,size_t,size_t,const GAOption*){};
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
            _otherOptFun=[](ArgsType*,std::list<Gene>*,size_t,size_t,const GAOption*){};
        } else {
            _otherOptFun=_ooF;
        }

        for(auto & i : _population) {
            _initializeFun(&i.self,&_args);
        }
        _eliteIt=_population.begin();
        _recording.clear();
        _recording.reserve(_option.maxGenerations+1);
    }


    ///start to solve
    virtual void run() {
        _generation=0;
        _failTimes=0;
        while(true) {
            _generation++;
            calculateAll();
            select();
            _recording.emplace_back(_eliteIt->fitness());
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
    typename std::list<Gene>::iterator eliteIt() const {
        return _eliteIt;
    }
    ///result
    const Var & result() const {
        return _eliteIt->self;
    }
    ///the whole population
    const std::list<Gene> & population() const {
        return _population;
    }
    const std::vector<double> & recording() const {
        return _recording;
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
    typedef typename std::list<Gene>::iterator GeneIt;
    GeneIt _eliteIt;
    std::list<Gene> _population;
    GAOption _option;

    size_t _generation;
    size_t _failTimes;
    std::vector<double> _recording;

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

        const double prevEliteFitness=_eliteIt->fitness();
        std::vector<GeneIt> iterators;
        iterators.clear();
        iterators.reserve(_population.size());
        auto GeneItCmp=[](GeneIt a,GeneIt b) {
            return isBetter(a->_Fitness,b->_Fitness);
        };

        for(auto it=_population.begin();it!=_population.end();++it) {
            iterators.emplace_back(it);
        }

        std::sort(iterators.begin(),iterators.end(),GeneItCmp);
        
        while(_population.size()>_option.populationSize) {
            _population.erase(iterators.back());
            iterators.pop_back();
        }

        GeneIt curBest=iterators.front();
        if(!isBetter(curBest->fitness(),prevEliteFitness)) {
            _failTimes++;
            _eliteIt=curBest;
        }
        else {
            _failTimes=0;
            _eliteIt=curBest;
        }

        _population.emplace_back(*_eliteIt);

    }

    virtual void crossover() {
        std::vector<GeneIt> crossoverQueue;
        crossoverQueue.clear();
        crossoverQueue.reserve(_population.size());

        for(GeneIt it=_population.begin();it!=_population.end();it++) {
            if(it==_eliteIt) {
                continue;
            }

            if(randD()<=_option.crossoverProb) {
                crossoverQueue.emplace_back(it);
            }
        }

        std::random_shuffle(crossoverQueue.begin(),crossoverQueue.end());

        if(crossoverQueue.size()%2==1) {
            crossoverQueue.pop_back();
        }

        while(crossoverQueue.empty()) {
            GeneIt a,b;
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
}
#endif // GENETIC_H
