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

#ifndef Heu_SOGA_H
#define Heu_SOGA_H
#include "./GABase.hpp"


namespace Heu
{

/**
   *  @brief Single object genetic algorithm solver(real fitness value).
   *
   *  @tparam Var_t  Type of decisition variable.
   *  @tparam isGreaterBetter Whether greater fitness value means better.
   *  @tparam Record  Whether the solver records fitness changelog.
   *  @tparam Args_t  Type of other parameters.
  */
template<typename Var_t,
    FitnessOption isGreaterBetter=FITNESS_LESS_BETTER,
    RecordOption Record=DONT_RECORD_FITNESS,
    class Args_t=void,
         typename GAAbstract<Var_t,double,Args_t>::initializeFun _iFun_=nullptr,
         typename GAAbstract<Var_t,double,Args_t>::fitnessFun _fFun_=nullptr,
         typename GAAbstract<Var_t,double,Args_t>::crossoverFun _cFun_=nullptr,
         typename GAAbstract<Var_t,double,Args_t>::mutateFun _mFun_=nullptr>
class SOGA : public GABase<Var_t,double,Record,Args_t,_iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    SOGA() {

    };
    using Base_t = GABase<Var_t,double,Record,Args_t,_iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_GABASE_TYPES

    virtual double bestFitness() const {
        return _eliteIt->_Fitness;
    }

    inline const Var_t & result() const {
        return _eliteIt->self;
    }

    void initializePop() {
        Base_t::initializePop();
        this->_eliteIt=this->_population.begin();
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
        iterators.reserve(this->_population.size());
        auto GeneItCmp=[](GeneIt_t a,GeneIt_t b) {
            return isBetter(a->_Fitness,b->_Fitness);
        };

        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            iterators.emplace_back(it);
        }

        std::sort(iterators.begin(),iterators.end(),GeneItCmp);
        
        while(this->_population.size()>this->_option.populationSize) {
            this->_population.erase(iterators.back());
            iterators.pop_back();
        }

        GeneIt_t curBest=iterators.front();
        if(!isBetter(curBest->_Fitness,prevEliteFitness)) {
            this->_failTimes++;
            _eliteIt=curBest;
        }
        else {
            this->_failTimes=0;
            _eliteIt=curBest;
        }

        this->_population.emplace_back(*_eliteIt);

    }

    virtual void mutate() {
        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            if(it==_eliteIt) {
                continue;
            }
            if(randD()<=this->_option.mutateProb) {

                if constexpr (std::is_same<Args_t,void>::value)
                        this->runmFun(&it->self);
                else
                    this->runmFun(&it->self,&this->_args);

                it->setUncalculated();
            }
        }
    }
};
}
#endif // SOGA_H
