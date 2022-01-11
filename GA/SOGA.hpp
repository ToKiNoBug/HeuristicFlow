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

#ifndef SOGA_H
#define SOGA_H
#include "./GABase.hpp"


namespace OptimT
{

///single object genetic algorithm solver(real fitness value)
template<typename Var_t,bool isGreaterBetter,bool Record,class ...Args>
class SOGA : public GABase<Var_t,double,Record,Args...>
{
public:
    SOGA() {
        //initialization for function ptrs
        Base_t::_initializeFun=
                [](Var_t*,const ArgsType*){};
        Base_t::_fitnessFun=
                [](const Var_t*,const ArgsType*,double*){};
        Base_t::_crossoverFun=
                [](const Var_t*,const Var_t*,Var_t*,Var_t*,const ArgsType*){};
        Base_t::_mutateFun=
                [](Var_t*,const ArgsType*){};
        Base_t::_otherOptFun=
                [](ArgsType*,std::list<typename Base_t::Gene>*,size_t,size_t,const GAOption*){};

  };
    using Base_t = GABase<Var_t,double,Record,Args...>;
    OPTIMT_MAKE_GABASE_TYPES

    virtual double bestFitness() const {
        return _eliteIt->_Fitness;
    }

    const Var_t & result() const {
        return _eliteIt->self;
    }

    virtual void initialize(
                            initializeFun _iFun,
                            fitnessFun _fFun,
                            crossoverFun _cFun,
                            mutateFun _mFun,
                            otherOptFun _ooF=nullptr,
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
            Base_t::_otherOptFun=[]
                    (ArgsType*,std::list<typename Base_t::Gene>*,size_t,size_t,const GAOption*){};
        } else {
            Base_t::_otherOptFun=_ooF;
        }

        for(auto & i : Base_t::_population) {
            Base_t::_initializeFun(&i.self,&Base_t::args());
            i.setUncalculated();
        }
        _eliteIt=Base_t::_population.begin();
    }

protected:

    GeneIt_t _eliteIt;

    static inline bool isBetter(double A,double B) {
        if(isGreaterBetter) {
            return A>B;
        }
        return A<B;
    }

    virtual void select() {

        const double prevEliteFitness=_eliteIt->_Fitness;
        std::vector<GeneIt_t> iterators;
        iterators.clear();
        iterators.reserve(Base_t::_population.size());
        auto GeneItCmp=[](GeneIt_t a,GeneIt_t b) {
            return isBetter(a->_Fitness,b->_Fitness);
        };

        for(auto it=Base_t::_population.begin();it!=Base_t::_population.end();++it) {
            iterators.emplace_back(it);
        }

        std::sort(iterators.begin(),iterators.end(),GeneItCmp);
        
        while(Base_t::_population.size()>Base_t::_option.populationSize) {
            Base_t::_population.erase(iterators.back());
            iterators.pop_back();
        }

        GeneIt_t curBest=iterators.front();
        if(!isBetter(curBest->_Fitness,prevEliteFitness)) {
            Base_t::_failTimes++;
            _eliteIt=curBest;
        }
        else {
            Base_t::_failTimes=0;
            _eliteIt=curBest;
        }

        Base_t::_population.emplace_back(*_eliteIt);

    }

    virtual void mutate() {
        for(auto it=Base_t::_population.begin();it!=Base_t::_population.end();++it) {
            if(it==_eliteIt) {
                continue;
            }
            if(OtGlobal::randD()<=Base_t::_option.mutateProb) {
                Base_t::_mutateFun(&it->self,&Base_t::args());
                it->setUncalculated();
            }
        }
    }
};
}
#endif // SOGA_H
