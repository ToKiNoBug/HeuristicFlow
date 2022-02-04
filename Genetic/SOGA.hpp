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

#ifndef OptimT_SOGA_H
#define OptimT_SOGA_H
#include "./GABase.hpp"


namespace OptimT
{

/**
   *  @brief Single object genetic algorithm solver(real fitness value).
   *
   *  @tparam Var_t  Type of decisition variable.
   *  @tparam isGreaterBetter Whether greater fitness value means better.
   *  @tparam Record  Whether the solver records fitness changelog.
   *  @tparam Args_t  Type of other parameters.
  */
template<typename Var_t,FitnessOption isGreaterBetter,RecordOption Record,class Args_t=void>
class SOGA : public GABase<Var_t,double,Record,Args_t>
{
public:
    SOGA() {

    };
    using Base_t = GABase<Var_t,double,Record,Args_t>;
    OptimT_MAKE_GABASE_TYPES

    virtual double bestFitness() const {
        return _eliteIt->_Fitness;
    }

    inline const Var_t & result() const {
        return _eliteIt->self;
    }

    virtual void customOptAfterInitialization() {
        Base_t::customOptAfterInitialization();
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
            if(randD()<=Base_t::_option.mutateProb) {
                Base_t::_mutateFun(&it->self,&Base_t::args());
                it->setUncalculated();
            }
        }
    }
};
}
#endif // SOGA_H
