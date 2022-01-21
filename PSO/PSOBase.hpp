#ifndef PSOBASE_HPP
#define PSOBASE_HPP

#include "PSOOption.hpp"

#ifndef OptimT_NO_OUTPUT
#include <iostream>
#endif

namespace OptimT {

template<class Var_t,class Fitness_t,RecordOption Record,class...Args>
class PSOBase
{
public:
    using Args_t = std::tuple<Args...>;
    using iFun_t = void(*)(Var_t * pos,Var_t * velocity,
        const Var_t * pMin,const Var_t * pMax,const Var_t * vMax,const Args_t *);
    using fFun_t = void(*)(const Var_t * ,const Args_t *,Fitness_t *);

    class Particle
    {
    public:
        Var_t position;
        Var_t velocity;
        bool isCalculated;
        Fitness_t fitness;
        void setUncalculated() {
            isCalculated=false;
        }
    };
    using ooFun_t = void(*)
        (Args_t*,std::vector<Particle>*,size_t generation,size_t failTimes,PSOOption*);

public:
    PSOBase() {};
    virtual ~PSOBase() {};

    const PSOOption & option() const {
        return _option;
    }

    const Args_t & args() const {
        return _args;
    }

    size_t generation() const {
        return _generation;
    }

    size_t failTimes() const {
        return _failTimes;
    }

    const Var_t & posMin() const {
        return _posMin;
    }

    const Var_t & posMax() const {
        return _posMax;
    }

    const Var_t & velocityMax() const {
        return _velocityMax;
    }

    iFun_t initializeFun() const {
        return _iFun;
    }

    fFun_t fitnessFun() const {
        return _fFun;
    }

    ooFun_t otherOperationFun() const {
        return _ooFun;
    }

    const std::vector<Particle> & population() const {
        return _population;
    }

    const Particle * populationBest() const {
        return pBest;
    }

    const Particle * globalBest() const {
        return &gBest;
    }

    void setPVRange(const Var_t & pMin,
            const Var_t & pMax,
            const Var_t & vMax) {
        _posMin=pMin;
        _posMax=pMax;
        _velocityMax=vMax;
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
            i.setUncalculated();
        }
        _generation=0;
        _failTimes=0;
    }

    virtual void run() {
        _generation=0;
        _failTimes=0;

        while(true) {
            _generation++;
            calculateAll();
            updatePGBest();
            if(_generation>_option.maxGeneration) {
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
                updatePopulation();
        }
        _generation--;

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

    const Particle * pBest;
    Particle gBest;

    virtual void calculateAll() {
#ifdef OptimT_DO_PARALLELIZE
        static const uint32_t thN=OtGlobal::threadNum();
        std::vector<Particle*> tasks;
        tasks.resize(0);
        tasks.reserve(_population.size());
        for(Particle & i : _population) {
            if(i._isCalculated) {
                continue;
            }
            tasks.emplace_back(&i);
        }
#pragma omp parallel for
        for(uint32_t begIdx=0;begIdx<thN;begIdx++) {
            for(uint32_t i=begIdx;i<tasks.size();i+=thN) {
                Particle * ptr=tasks[i];
                _fitnessFun(&ptr->position,&_args,&ptr->_Fitness);
                ptr->isCalculated=true;
            }
        }
#else
        for(Particle & i : _population) {
            if(i.isCalculated) {
                continue;
            }
            _Fun(&i.position,&_args,&i._Fitness);
            i.isCalculated=true;
        }
#endif
    }

    virtual void updatePGBest()=0;

    virtual void updatePopulation()=0;

};



template<class Var_t,class Fitness_t,class...Args>
class PSOBase<Var_t,Fitness_t,RECORD_FITNESS,Args...>
        : public PSOBase<Var_t,Fitness_t,DONT_RECORD_FITNESS,Args...>
{

};

}
#endif // PSOBASE_HPP
