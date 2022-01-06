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
#include "./GABase.h"

#include <vector>
#include <tuple>
#include <algorithm>

namespace OptimT
{

///single object genetic algorithm solver(real fitness value)
template<typename Var_t,bool isGreaterBetter,bool Record,class ...Args>
class SOGA : public GABase<Var_t,double,Record,Args...>
{
public:
    SOGA() {
        //initialization for function ptrs has been finished in base class constructor
  };
    typedef GABase<Var_t,double,Record,Args...> Base_t;
protected:

    static inline bool isBetter(double A,double B) {
        if(isGreaterBetter) {
            return A>B;
        }
        return A<B;
    }

    virtual void select() {

        const double prevEliteFitness=Base_t::_eliteIt->_Fitness;
        typedef typeof(Base_t::_population.begin()) GeneIt_t;
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
            Base_t::_eliteIt=curBest;
        }
        else {
            Base_t::_failTimes=0;
            Base_t::_eliteIt=curBest;
        }

        Base_t::_population.emplace_back(*Base_t::_eliteIt);

    }
};
}
#endif // SOGA_H
