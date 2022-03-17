// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
            if
        #if __cplusplus >=201703L
                constexpr
        #endif
        (fOpt==FITNESS_GREATER_BETTER) {
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
        if
        #if __cplusplus >=201703L
                constexpr
        #endif
        (fOpt==FITNESS_GREATER_BETTER) {
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
