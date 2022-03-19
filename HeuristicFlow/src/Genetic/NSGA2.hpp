// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef Heu_NSGA2BASE_HPP
#define Heu_NSGA2BASE_HPP
#include "NSGABase.hpp"

namespace Heu
{


/**
 * @brief NSGA2 MOGA solver. Suitable for not too many objectives.
 * 
 * @tparam Var_t Type of decisition variable.
 * @tparam ObjNum Numbers of objectives.
 * @tparam isGreaterBetter Whether greater fitness value means better.
 * @tparam Record Whether the solver records fitness changelog.
 * @tparam Args_t Type of other parameters.
 * @tparam _iFun_ Compile-time iFun, use nullptr for runtime
 * @tparam _fFun_ Compile-time fFun, use nullptr for runtime
 * @tparam _cFun_ Compile-time cFun, use nullptr for runtime
 * @tparam _mFun_ Compile-time mFun, use nullptr for runtime
 */
template<typename Var_t,
         size_t ObjNum,
         FitnessOption fOpt=FITNESS_LESS_BETTER,
         RecordOption rOpt=DONT_RECORD_FITNESS,
         class Args_t=void,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::initializeFun _iFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::fitnessFun _fFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::crossoverFun _cFun_=nullptr,
         typename internal::GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::mutateFun _mFun_=nullptr>
class NSGA2
    :public internal::NSGABase<Var_t,ObjNum,fOpt,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>
{
    using Base_t = internal::NSGABase<Var_t,ObjNum,fOpt,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
private:
public:
    NSGA2() {};
    virtual ~NSGA2() {};
    Heu_MAKE_NSGABASE_TYPES

    /** @brief temporary struct to store infos when selection
     */
    struct infoUnit2 : public infoUnitBase_t
    {
    public:
        /** @brief whether this gene is selected
        */
        double congestion;
    };

protected:
    static bool compareByCongestion(const infoUnitBase_t * A,const infoUnitBase_t * B) {
        if(A==B)
            return false;
        return (static_cast<const infoUnit2*>(A)->congestion) > (static_cast<const infoUnit2*>(B)->congestion);
    }

    template<int64_t objIdx>
    static bool compareByFitness(const infoUnitBase_t * A,const infoUnitBase_t * B) {
#ifndef Heu_NO_STATICASSERT
    static_assert(objIdx>=0,"Heu : Invalid comparison flag");
#endif
        if(A==B) 
            return false;
        ///compare by fitness on single objective
        return A->iterator->_Fitness[objIdx]<B->iterator->_Fitness[objIdx];
    }

    ///fast nondominated sorting
    virtual void select() {
        using cmpFun_t = bool(*)(const infoUnitBase_t * ,const infoUnitBase_t * );
        static const size_t objCapacity=
                (ObjNum==0)?(Heu_MOGA_MaxRunTimeObjNum):ObjNum;
        static const std::array<cmpFun_t,objCapacity> fitnessCmpFuns
                =expand<0,objCapacity-1>();

        const size_t popSizeBefore=this->_population.size();
        std::vector<infoUnit2> pop;
        pop.clear();pop.reserve(popSizeBefore);
        this->sortSpace.reserve(2*this->_option.populationSize);
        this->sortSpace.resize(popSizeBefore);

        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            pop.emplace_back();
            pop.back().iterator=it;
            pop.back().congestion=0;
        }

        //make sortspace
        for(size_t i=0;i<popSizeBefore;i++) {
            this->sortSpace[i]=pop.data()+i;
        }
        
        this->calculateDominatedNum();

        this->divideLayers();
        
        const size_t PFSize=this->pfLayers.front().size();
        if(PFSize<=this->_option.populationSize)
            this->updatePF((const infoUnitBase_t **)(this->pfLayers.front().data()),
                         this->pfLayers.front().size());


        std::unordered_set<infoUnit2 *> selected;
        selected.reserve(this->_option.populationSize);
        bool needCongestion=true;
        while(true) {
            //don't need to calculate congestion
            if(selected.size()==this->_option.populationSize) {
                needCongestion=false;
                break;
            }
            //need to calculate congestion
            if(selected.size()+this->pfLayers.front().size()>this->_option.populationSize) {
                needCongestion=true;
                break;
            }
            //emplace every element of this layer into selected
            for(const auto i : this->pfLayers.front()) {
                selected.emplace(static_cast<infoUnit2 *>(i));
            }
            this->pfLayers.pop_front();
        }
        //calculate congestion
        if(needCongestion) {
            for(size_t objIdx=0;objIdx<this->objectiveNum();objIdx++) {

                std::sort((infoUnit2 **)(this->sortSpace.data()),
                            (infoUnit2 **)(this->sortSpace.data()+this->sortSpace.size()),
                          fitnessCmpFuns[objIdx]);

                const double scale=std::abs(this->sortSpace.front()->iterator->_Fitness[objIdx]
                        -this->sortSpace.back()->iterator->_Fitness[objIdx])
                    +1e-10;

                static_cast<infoUnit2 *>(this->sortSpace.front())->congestion=pinfD;
                static_cast<infoUnit2 *>(this->sortSpace.back())->congestion=pinfD;

                //calculate congestion on single object
                for(size_t idx=1;idx<popSizeBefore-1;idx++) {

                    static_cast<infoUnit2 *>(this->sortSpace[idx])->congestion
                            +=std::abs(
                            this->sortSpace[idx-1]->iterator->_Fitness[objIdx]
                            -this->sortSpace[idx+1]->iterator->_Fitness[objIdx]
                            )/scale;
                }
            } // end sort on objIdx
            
            //sort by congestion in the undetermined set
            std::sort((infoUnit2 **)(this->pfLayers.front().data()),
                      (infoUnit2 **)(this->pfLayers.front().data()+this->pfLayers.front().size()),
                      compareByCongestion);

            size_t idx=0;
            while(selected.size()<this->_option.populationSize) {
                selected.emplace(static_cast<infoUnit2 *>(this->pfLayers.front()[idx]));
                idx++;
            }

            

        } // end applying congestion


        
        //erase unselected
        for(auto & i : pop) {
            if(selected.find(&i)==selected.end()) {
                this->_population.erase(i.iterator);
            }
        }

        if(PFSize>this->_option.populationSize) {
            std::vector<const infoUnitBase_t*> PF;
            PF.reserve(selected.size());
            for(auto i : selected) {
                PF.emplace_back(i);
            }
            this->updatePF(PF.data(),PF.size());
        }


        this->pfLayers.clear();
        this->sortSpace.clear();
        
    }


private:
    //some template metaprogramming to make a function pointer array as below:
    //universialCompareFun<0>,universialCompareFun<1>,...,universialCompareFun<ObjNum-1>
    using fun_t = bool(*)(const infoUnitBase_t * ,const infoUnitBase_t * );
    
    template<int64_t cur,int64_t max>
    struct expandStruct
    {
        static void expand(fun_t * dst) {
            *dst=compareByFitness<cur>;
            expandStruct<cur+1,max>::expand(dst+1);
        }
    };

    template<int64_t max>
    struct expandStruct<max,max>
    {
        static void expand(fun_t * dst) {
            *dst=compareByFitness<max>;
        }
    };

    template<int64_t beg,int64_t end>
    std::array<fun_t,end-beg+1> expand() {
        std::array<fun_t,end-beg+1> funs;
        expandStruct<beg,end>::expand(funs.data());
        return funs;
    }


};


}   // Heu

#endif //   NSGA2BASE_HPP
