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

#ifndef OptimT_NSGA3_HPP
#define OptimT_NSGA3_HPP

#include "NSGA3Base.hpp"

namespace OptimT {

template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        RecordOption rOpt,
        PFOption pfOpt,
        ReferencePointOption rpOpt,
        class Args_t>
class NSGA3 : public NSGA3Base<Var_t,ObjNum,DVO,rOpt,pfOpt,rpOpt,Args_t>
{
public:
    NSGA3() {};
    virtual ~NSGA3() {};
    using Base_t = NSGA3Base<Var_t,ObjNum,DVO,rOpt,pfOpt,rpOpt,Args_t>;
    OptimT_MAKE_NSGA3ABSTRACT_TYPES

    virtual Fitness_t bestFitness() const {
        Fitness_t best=this->_population.front()._Fitness;
        for(const Gene * i : this->_pfGenes) {
            for(size_t objIdx=0;objIdx<i->_Fitness.size();objIdx++) {
                best[objIdx]=std::min(best[objIdx],i->_Fitness[objIdx]);
            }
        }
        return best;
    }

    
    void initializePop() {
        this->makeReferencePoses();
        Base_t::initializePop();
    }
};

}   //  namespace OptimT

#endif  //  OptimT_NSGA3_HPP
