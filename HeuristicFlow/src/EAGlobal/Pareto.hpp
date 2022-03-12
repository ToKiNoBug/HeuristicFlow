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

#ifndef Heu_PARETO_HPP
#define Heu_PARETO_HPP

#include "../Global/Enumerations.hpp"
#include "../Global/Types.hpp"
#include "../Global/Constants.hpp"

namespace Heu {
template<size_t ObjNum,DoubleVectorOption dvo,FitnessOption fOpt>
struct Pareto
{
    using Fitness_t = FitnessVec_t<dvo,ObjNum>;
    static bool isStrongDominate(const Fitness_t * A,const Fitness_t * B) {
        if(A==B) return false;
        uint32_t notWorseNum=0,betterNum=0;
        for(size_t objIdx=0;objIdx<A->size();objIdx++) {
            if constexpr (fOpt==FITNESS_GREATER_BETTER) {
                notWorseNum+=((*A)[objIdx]>=(*B)[objIdx]);
                betterNum+=((*A)[objIdx]>(*B)[objIdx]);
            }
            else {
                notWorseNum+=((*A)[objIdx]<=(*B)[objIdx]);
                betterNum+=((*A)[objIdx]<(*B)[objIdx]);
            }
        }
        if(notWorseNum<A->size())
            return false;
        return betterNum>0;
    }

};

#ifdef EIGEN_CORE_H

template<size_t ObjNum,FitnessOption fOpt>
struct Pareto<ObjNum,DoubleVectorOption::Eigen,fOpt>
{
    using Fitness_t = EigenVecD_t<ObjNum>;
    static bool isStrongDominate(const Fitness_t * A,const Fitness_t * B) {
        bool isNotWorse,isBetter;
        if constexpr (fOpt==FITNESS_GREATER_BETTER) {
            isNotWorse=((*A)>=(*B)).all();
            isBetter=((*A)>(*B)).any();
        }
        else {
            isNotWorse=((*A)<=(*B)).all();
            isBetter=((*A)<(*B)).any();
        }
        return isNotWorse&&isBetter;
    }
};

#endif  //  EIGEN_CORE_H

}   //  namespace Heu


#endif // Heu_PARETO_HPP
