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

#ifndef Heu_NSGA2BASE_HPP
#define Heu_NSGA2BASE_HPP

#include "NSGABase.hpp"

namespace Heu {

enum CompareOption : int64_t {
    CompareByCongestion=-1,
    CompareByDominantedBy=-2
};

/**
 * @brief The NSGA2Base class implements most NSGA-II functions
 * and can be inherited inorder to boost with Eigen
 */
template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class Args_t,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::mutateFun _mFun_>
class NSGA2Base
    :public NSGABase<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    NSGA2Base() {
        _ccFun=default_ccFun_liner;
    };
    virtual ~NSGA2Base() {};

    using Base_t = NSGABase<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_NSGABASE_TYPES

    using congestComposeFun = double(*)(const Fitness_t *);

    inline void setCongestComposeFun(congestComposeFun __ccFun=default_ccFun_liner) {
            _ccFun=__ccFun;
    }

    static double default_ccFun_liner(const Fitness_t * f) {
        double result=f->operator[](0);
        for(size_t objIdx=1;objIdx<((ObjNum==Dynamic)?f->size():ObjNum);objIdx++) {
            result+=f->operator[](objIdx);
        }
        return result;
    };

    static double default_ccFun_sphere(const Fitness_t * f) {
        double result=OT_square(f->operator[](0));
        for(size_t objIdx=1;objIdx<((ObjNum==Dynamic)?f->size():ObjNum);objIdx++) {
            result+=OT_square(f->operator[](objIdx));
        }
        return std::sqrt(result);
    }

    static double default_ccFun_max(const Fitness_t * f) {
        double result=f->operator[](0);
        for(size_t objIdx=1;objIdx<((ObjNum==Dynamic)?f->size():ObjNum);objIdx++) {
            result=std::max(f->operator[](objIdx),result);
        }
        return result;
    }

    /**
     * @brief compose congestion using p-th minkowski distance
     * 
     * @tparam p power
     * @param f parital congestion value
     * @return double congestion value
     */
    template<int64_t p>
    static double default_ccFun_powered(const Fitness_t * f) {
        double result=0;
        for(size_t objIdx=0;objIdx<((ObjNum==Dynamic)?f->size():ObjNum);objIdx++) {
            result+=power<p>(f->operator[](objIdx));
        }
        return std::pow(result,1.0/p);
    }

    inline void initializePop() {
        if(_ccFun==nullptr) {
            setCongestComposeFun();
        }
        Base_t::initializePop();
    }

    /**
     * @brief calculate ideal point
     * 
     * @return Fitness_t ideal point
     */
    virtual Fitness_t bestFitness() const {
        Fitness_t best=this->_population.front()._Fitness;
        for(const Gene & i : this->_population) {
            if(fOpt) {
                for(size_t objIdx=0;objIdx<this->objectiveNum();objIdx++) {
                    best[objIdx]=std::max(best[objIdx],i._Fitness[objIdx]);
                }
            } else {
                for(size_t objIdx=0;objIdx<this->objectiveNum();objIdx++) {
                    best[objIdx]=std::min(best[objIdx],i._Fitness[objIdx]);
                }
            }
        }
        return best;
    }

    /** @brief temporary struct to store infos when selection
     */
    struct infoUnit2 : public infoUnitBase_t
    {
    public:
        /** @brief whether this gene is selected
        */
        Fitness_t congestion;
    };

protected:
    congestComposeFun _ccFun;

    template<int64_t objIdx>
    static bool universialCompareFun(const infoUnit2 * A,const infoUnit2 * B) {
#ifndef Heu_NO_STATICASSERT
    static_assert(std::integral_constant<bool,
        ((objIdx>=0)
        ||(objIdx==CompareOption::CompareByCongestion)
        ||(objIdx==CompareOption::CompareByDominantedBy))>::value,
    "Heu : Invalid compare flag");
#endif
        if(A==B) return false;

#if __cplusplus >= 201703L
        ///compare by congestion
        if constexpr (objIdx==CompareByCongestion) {
            return A->congestion[0]>B->congestion[0];
        }
        ///compare by this->pfLayers
        if constexpr(objIdx==CompareByDominantedBy) {
            return A->domainedByNum<B->domainedByNum;
        }
#else   //  these if constexpr will be done by compiler's optimization
        ///compare by congestion
        if (objIdx==CompareByCongestion) {
            return A->congestion[0]>B->congestion[0];
        }
        ///compare by this->pfLayers
        if (objIdx==CompareByDominantedBy) {
            return A->domainedByNum<B->domainedByNum;
        }
#endif  //  #if __cplusplus >= 201703L

        ///compare by fitness on single objective
        return A->iterator->_Fitness[objIdx]<B->iterator->_Fitness[objIdx];
    }

    ///fast nondominated sorting
    virtual void select() {
        using cmpFun_t = bool(*)(const infoUnit2 * ,const infoUnit2 * );
        static const size_t objCapacity=
                (ObjNum==0)?(Heu_MOGA_MaxRunTimeObjNum):ObjNum;
        static const std::array<cmpFun_t,objCapacity> fitnessCmpFuns
                =expand<0,objCapacity-1>();

        const size_t popSizeBefore=this->_population.size();
        std::vector<infoUnit2> pop;
        pop.clear();pop.reserve(popSizeBefore);

        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            pop.emplace_back();
            pop.back().iterator=it;

            initializeSize<ObjNum>::
                    template resize<Fitness_t>(&pop.back().congestion,this->objectiveNum());

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
                selected.emplace((infoUnit2 *)i);
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
                    +1e-100;

                ((infoUnit2 *)this->sortSpace.front())->congestion[objIdx]=Heu::pinfD;
                ((infoUnit2 *)this->sortSpace.back())->congestion[objIdx]=Heu::pinfD;

                //calculate congestion on single object
                for(size_t idx=1;idx<popSizeBefore-1;idx++) {

                    ((infoUnit2 *)this->sortSpace[idx])->congestion[objIdx]=std::abs(
                                this->sortSpace[idx-1]->iterator->_Fitness[objIdx]
                               -this->sortSpace[idx+1]->iterator->_Fitness[objIdx]
                                )/scale;
                }
            } // end sort on objIdx

            for(infoUnit2 & i : pop) {
                //store final congestion at the first congestion
                i.congestion[0]=_ccFun(&i.congestion);
            }

            std::sort((infoUnit2 **)(this->sortSpace.data()),
                      (infoUnit2 **)(this->sortSpace.data()+this->sortSpace.size()),
                      universialCompareFun<CompareByCongestion>);
            size_t idx=0;
            while(selected.size()<this->_option.populationSize) {
                selected.emplace((infoUnit2 *)this->pfLayers.front()[idx]);
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
    }


private:
    //some template metaprogramming to make a function pointer array as below:
    //universialCompareFun<0>,universialCompareFun<1>,...,universialCompareFun<ObjNum-1>
    using fun_t = bool(*)(const infoUnit2 * ,const infoUnit2 * );
    template<int64_t cur,int64_t max>
    struct expandStruct
    {
        static void expand(fun_t * dst) {
            *dst=universialCompareFun<cur>;
            expandStruct<cur+1,max>::expand(dst+1);
        }
    };

    template<int64_t max>
    struct expandStruct<max,max>
    {
        static void expand(fun_t * dst) {
            *dst=universialCompareFun<max>;
        }
    };

    template<int64_t beg,int64_t end>
    std::array<fun_t,end-beg+1> expand() {
        std::array<fun_t,end-beg+1> funs;
        expandStruct<beg,end>::expand(funs.data());
        return funs;
    }


};

}

#endif //   NSGA2BASE_HPP
