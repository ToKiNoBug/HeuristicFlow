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

#ifndef OptimT_NSGA2BASE_HPP
#define OptimT_NSGA2BASE_HPP

#include "NSGABase.hpp"

namespace OptimT {

#ifdef OptimT_NSGA2_DO_PARALLELIZE
#ifndef OptimT_DO_PARALLELIZE
#error You allowed parallelize in NSGA2 but not on global.  \
    Macro OptimT_NSGA2_DO_PARALLELIZE can only be defined when OptimT_DO_PARALLELIZE is defined.
#endif
#endif


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
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class ...Args>
class NSGA2Base
    :public NSGABase<Var_t,ObjNum,Fitness_t,fOpt,rOpt,pfOpt,Args...>
{
public:
    NSGA2Base() {
        _ccFun=default_ccFun_liner;
    };
    virtual ~NSGA2Base() {};

    using Base_t = NSGABase<Var_t,ObjNum,Fitness_t,fOpt,rOpt,pfOpt,Args...>;
    OptimT_MAKE_NSGABASE_TYPES

    using congestComposeFun = double(*)(const Fitness_t *,const ArgsType*);

    inline void setCongestComposeFun(congestComposeFun __ccFun=default_ccFun_liner) {
            _ccFun=__ccFun;
    }

    static double default_ccFun_liner(const Fitness_t * f,const ArgsType*) {
        double result=f->operator[](0);
        for(size_t objIdx=1;objIdx<((ObjNum==Dynamic)?f->size():ObjNum);objIdx++) {
            result+=f->operator[](objIdx);
        }
        return result;
    };

    static double default_ccFun_sphere(const Fitness_t * f,const ArgsType*) {
        double result=OT_square(f->operator[](0));
        for(size_t objIdx=1;objIdx<((ObjNum==Dynamic)?f->size():ObjNum);objIdx++) {
            result+=OT_square(f->operator[](objIdx));
        }
        return std::sqrt(result);
    }

    static double default_ccFun_max(const Fitness_t * f,const ArgsType*) {
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
    static double default_ccFun_powered(const Fitness_t * f,const ArgsType*) {
        double result=0;
        for(size_t objIdx=0;objIdx<((ObjNum==Dynamic)?f->size():ObjNum);objIdx++) {
            result+=power<p>(f->operator[](objIdx));
        }
        return std::pow(result,1.0/p);
    }
    /**
     * @brief calculate ideal point
     * 
     * @return Fitness_t ideal point
     */
    virtual Fitness_t bestFitness() const {
        Fitness_t best=Base_t::_population.front()._Fitness;
        for(const Gene & i : Base_t::_population) {
            if(fOpt) {
                for(size_t objIdx=0;objIdx<Base_t::objectiveNum();objIdx++) {
                    best[objIdx]=std::max(best[objIdx],i._Fitness[objIdx]);
                }
            } else {
                for(size_t objIdx=0;objIdx<Base_t::objectiveNum();objIdx++) {
                    best[objIdx]=std::min(best[objIdx],i._Fitness[objIdx]);
                }
            }
        }
        return best;
    }

    /** @brief temporary struct to store infos when selection
     */
    struct infoUnit:public Base_t::infoUnitBase
    {
    public:
        /** @brief whether this gene is selected
        */
        bool isSelected;
        Fitness_t congestion;
    };

protected:
    congestComposeFun _ccFun;

    virtual void customOptWhenInitialization() {
        this->prevFrontSize=-1;
        this->_pfGenes.clear();
        this->_pfGenes.reserve(Base_t::_option.populationSize*2);
        if(_ccFun==nullptr) {
            setCongestComposeFun();
        }
    }

    template<int64_t objIdx>
    static bool universialCompareFun(const infoUnit * A,const infoUnit * B) {
#ifndef OptimT_NO_STATICASSERT
    static_assert(std::integral_constant<bool,
        ((objIdx>=0)
        ||(objIdx==CompareOption::CompareByCongestion)
        ||(objIdx==CompareOption::CompareByDominantedBy))>::value,
    "OptimTemplates : Invalid compare flag");
#endif
        if(A==B) return false;
        ///compare by congestion
        if(objIdx==CompareByCongestion) {
            return A->congestion[0]>B->congestion[0];
        }
        ///compare by paretoLayers
        if(objIdx==CompareByDominantedBy) {
            return A->domainedByNum<B->domainedByNum;
        }
        ///compare by fitness on single objective
        return A->iterator->_Fitness[objIdx]<B->iterator->_Fitness[objIdx];
    }

    ///fast nondominated sorting
    virtual void select() {
        using cmpFun_t = bool(*)(const infoUnit * ,const infoUnit * );
        static const size_t objCapacity=
                (ObjNum==0)?(OptimT_MOGA_RTObjNum_MaxObjNum):ObjNum;
        static const std::array<cmpFun_t,objCapacity> fitnessCmpFuns
                =expand<0,objCapacity-1>();

        const size_t popSizeBefore=Base_t::_population.size();
        std::vector<infoUnit> pop;
        pop.clear();pop.reserve(popSizeBefore);

        for(auto it=Base_t::_population.begin();it!=Base_t::_population.end();++it) {
            pop.emplace_back();
            pop.back().isSelected=false;
            pop.back().iterator=it;
        }

        std::vector<infoUnit*> sortSpace(popSizeBefore);
        //make sortspace
        for(size_t i=0;i<popSizeBefore;i++) {
            sortSpace[i]=pop.data()+i;
        }

        Base_t::calculateDominatedNum((infoUnitBase_t **)sortSpace.data(),popSizeBefore);

        //sort by paretoLayers
        std::sort(sortSpace.begin(),sortSpace.end(),
                universialCompareFun<CompareByDominantedBy>);

        std::list<std::vector<infoUnit *>> paretoLayers;
        //seperate them into layers
        {
            size_t unLayeredNum=popSizeBefore;
            size_t prevDomainedByNum=-1;
            for(const auto i : sortSpace) {
                if(i->domainedByNum!=prevDomainedByNum) {
                    paretoLayers.emplace_back();
                    paretoLayers.back().clear();
                    paretoLayers.back().reserve(unLayeredNum);
                    prevDomainedByNum=i->domainedByNum;
                }
                paretoLayers.back().emplace_back(i);
                unLayeredNum--;
            }
        }

        this->updatePF((const infoUnitBase_t **)paretoLayers.front().data(),
                         paretoLayers.front().size());


        std::queue<infoUnit *> selected;
        bool needCongestion=true;
        while(true) {
            //don't need to calculate congestion
            if(selected.size()==Base_t::_option.populationSize) {
                needCongestion=false;
                break;
            }
            //need to calculate congestion
            if(selected.size()+paretoLayers.front().size()>Base_t::_option.populationSize) {
                needCongestion=true;
                break;
            }
            //emplace every element of this layer into selected
            for(const auto i : paretoLayers.front()) {
                selected.emplace(i);
            }
            paretoLayers.pop_front();
        }

        //calculate congestion
        if(needCongestion) {
#ifdef OptimT_NSGA2_DO_PARALLELIZE
#pragma omp parallel for
#endif
            for(size_t objIdx=0;objIdx<Base_t::objectiveNum();objIdx++) {
                //if don't parallelize,  cursortSpace is only a reference to sortSpace;
                //otherwise it's copied to enable sorting concurrently
                std::vector<infoUnit*>
#ifndef OptimT_DO_PARALLELIZE
                        &
#endif
                        cursortSpace=sortSpace;

                std::sort(cursortSpace.begin(),cursortSpace.end(),fitnessCmpFuns[objIdx]);

                cursortSpace.front()->congestion[objIdx]=OptimT::pinfD;
                cursortSpace.back()->congestion[objIdx]=OptimT::pinfD;

                //calculate congestion on single object
                for(size_t idx=1;idx<popSizeBefore-1;idx++) {

                    cursortSpace[idx]->congestion[objIdx]=std::abs(
                                cursortSpace[idx-1]->iterator->_Fitness[objIdx]
                               -cursortSpace[idx+1]->iterator->_Fitness[objIdx]
                                );
                }
            } // end sort on objIdx

            for(infoUnit & i : pop) {
                //store final congestion at the first congestion
                i.congestion[0]=_ccFun(&i.congestion,&Base_t::args());
            }

            std::sort(paretoLayers.front().begin(),paretoLayers.front().end(),
                      universialCompareFun<CompareByCongestion>);
            size_t idx=0;
            while(selected.size()<Base_t::_option.populationSize) {
                selected.emplace(paretoLayers.front()[idx]);
                idx++;
            }

        } // end applying congestion

        //mark selected genes
        while(!selected.empty()) {
            selected.front()->isSelected=true;
            selected.pop();
        }
        //erase unselected
        for(infoUnit & i : pop) {
            if(!i.isSelected) {
                Base_t::_population.erase(i.iterator);
            }
        }
    }


private:
    //some template metaprogramming to make a function pointer array as below:
    //universialCompareFun<0>,universialCompareFun<1>,...,universialCompareFun<ObjNum-1>
    using fun_t = bool(*)(const infoUnit * ,const infoUnit * );
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
