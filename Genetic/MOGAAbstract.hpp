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

#ifndef Heu_MOGAABSTRACT_HPP
#define Heu_MOGAABSTRACT_HPP

#include "./GABase.hpp"
#include <queue>
#include <unordered_set>

namespace Heu {
    
///whether to protect pareto front when mutation or not
enum PFOption : unsigned char {
    PARETO_FRONT_DONT_MUTATE=true,
    PARETO_FRONT_CAN_MUTATE=false
};

/**
   *  @brief Base class for multi-objective genetic algorithm solver.
   *
   *  @tparam Var_t  Type of decisition variable.
   *  @tparam ObjNum Numbers of objectives.
   *  @tparam Fitness_t Type of fitness value.
   *  @tparam fOpt Whether greater fitness value means better.
   *  @tparam rOpt Whether the solver records fitness changelog.
   *  @tparam pfOpt Whether to protect the Pareto front from mutation.
   *  @tparam Args_t Type of other parameters.
  */
template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class Args_t>
class MOGAAbstract
    : public GABase<Var_t,FitnessVec_t<DVO,ObjNum>,rOpt,Args_t>
{
public:
    MOGAAbstract()  {};
    virtual ~MOGAAbstract() {};

    using Base_t = GABase<Var_t,FitnessVec_t<DVO,ObjNum>,rOpt,Args_t>;
    Heu_MAKE_GABASE_TYPES
    using Fitness_t = FitnessVec_t<DVO,ObjNum>;

    ///get pareto front in vec
    inline void paretoFront(std::vector<Fitness_t> & front) const {
        front.clear();  front.reserve(_pfGenes.size());
        for(const Gene* i : _pfGenes) {
            front.emplace_back(i->_Fitness);
        }
        return;
    }

    inline void paretoFront(std::vector<std::pair<const Var_t*,const Fitness_t*>> & front) const {
        front.clear();
        front.reserve(_pfGenes.size());
        for(const Gene* i : _pfGenes) {
            front.emplace_back(std::make_pair(&(i->self),&(i->_Fitness)));
        }
    }

    inline const std::unordered_set<const Gene*> & pfGenes() const {
        return _pfGenes;
    }

    inline void initializePop() {
        this->prevFrontSize=-1;
        this->_pfGenes.clear();
        this->_pfGenes.reserve(this->_option.populationSize*2);
        Base_t::initializePop();
    }

protected:
    size_t prevFrontSize;
    size_t prevPFCheckSum;
    std::unordered_set<const Gene*> _pfGenes;

    ///whether A strong domainates B
    static bool isStrongDomain(const Fitness_t * A,const Fitness_t * B) {
        if(A==B) return false;
        if constexpr (DVO!=Eigen) {
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
        else {
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
    } //isStrongDomain

    virtual size_t makePFCheckSum() const {
        std::vector<const Gene*> pfvec;
        pfvec.reserve(_pfGenes.size());
        for(auto i : _pfGenes) {
            pfvec.emplace_back(i);
        }
        std::sort(pfvec.begin(),pfvec.end());

        static const auto hashFun=_pfGenes.hash_function();
        size_t checkSum=hashFun(pfvec.front());
        for(size_t i=1;i<pfvec.size();i++) {
            checkSum^=hashFun(pfvec[i]);
        }
        return checkSum;
    }
    

    ///mutate operation
    virtual void mutate() {
        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            if(randD()<=this->_option.mutateProb) {
                if(pfOpt){
                    if(_pfGenes.find(&*it)!=_pfGenes.end()) {
                        continue;
                    }
                }

                Base_t::doMutate(this->_mFun,&it->self);

                it->setUncalculated();
            }
        }
    }

private:
#ifndef Heu_NO_STATICASSERT
    static_assert(std::integral_constant<bool,(ObjNum!=1)>::value,
    "HeuristicFlow : You assigned single objective in MOGA");

#ifndef EIGEN_CORE_H
    static_assert(DVO!=DoubleVectorOption::Eigen,
        "Include Eigen before using Eigen arrays as Fitness types");
#endif  //  EIGEN_CORE_H

#endif

};  // MOGAAbstract

}   //  Heu

#endif //   MOGAABSTRACT_HPP
