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
        FitnessOption fOpt=FITNESS_LESS_BETTER,
        RecordOption rOpt=DONT_RECORD_FITNESS,
        PFOption pfOpt=PARETO_FRONT_CAN_MUTATE,
        class Args_t=void>
class NSGA2Base
    :public NSGABase<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,Args_t>
{
public:
    NSGA2Base() {
        _ccFun=default_ccFun_liner;
    };
    virtual ~NSGA2Base() {};

    using Base_t = NSGABase<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,Args_t>;
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

    template<int64_t objIdx>
    static bool universialCompareFun(const infoUnit * A,const infoUnit * B) {
#ifndef Heu_NO_STATICASSERT
    static_assert(std::integral_constant<bool,
        ((objIdx>=0)
        ||(objIdx==CompareOption::CompareByCongestion)
        ||(objIdx==CompareOption::CompareByDominantedBy))>::value,
    "Heu : Invalid compare flag");
#endif
        if(A==B) return false;
        ///compare by congestion
        if constexpr(objIdx==CompareByCongestion) {
            return A->congestion[0]>B->congestion[0];
        }
        ///compare by paretoLayers
        if constexpr(objIdx==CompareByDominantedBy) {
            return A->domainedByNum<B->domainedByNum;
        }
        ///compare by fitness on single objective
        return A->iterator->_Fitness[objIdx]<B->iterator->_Fitness[objIdx];
    }

    ///fast nondominated sorting
    virtual void select() {
        using cmpFun_t = bool(*)(const infoUnit * ,const infoUnit * );
        static const size_t objCapacity=
                (ObjNum==0)?(Heu_MOGA_MaxRunTimeObjNum):ObjNum;
        static const std::array<cmpFun_t,objCapacity> fitnessCmpFuns
                =expand<0,objCapacity-1>();

        const size_t popSizeBefore=this->_population.size();
        std::vector<infoUnit> pop;
        pop.clear();pop.reserve(popSizeBefore);

        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            pop.emplace_back();
            pop.back().isSelected=false;
            pop.back().iterator=it;
            if constexpr (ObjNum==Dynamic) {
                if constexpr (DVO==DoubleVectorOption::Eigen) {
                    pop.back().congestion.resize(this->objectiveNum(),1);
                }
                else {
                    pop.back().congestion.resize(this->objectiveNum());
                }
            }
        }

        std::vector<infoUnit*> sortSpace(popSizeBefore);
        //make sortspace
        for(size_t i=0;i<popSizeBefore;i++) {
            sortSpace[i]=pop.data()+i;
        }

        this->calculateDominatedNum((infoUnitBase_t **)sortSpace.data(),popSizeBefore);

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
            if(selected.size()==this->_option.populationSize) {
                needCongestion=false;
                break;
            }
            //need to calculate congestion
            if(selected.size()+paretoLayers.front().size()>this->_option.populationSize) {
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
            for(size_t objIdx=0;objIdx<this->objectiveNum();objIdx++) {
                std::vector<infoUnit*> & cursortSpace=sortSpace;

                std::sort(cursortSpace.begin(),cursortSpace.end(),fitnessCmpFuns[objIdx]);

                cursortSpace.front()->congestion[objIdx]=Heu::pinfD;
                cursortSpace.back()->congestion[objIdx]=Heu::pinfD;

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
                i.congestion[0]=_ccFun(&i.congestion);
            }

            std::sort(paretoLayers.front().begin(),paretoLayers.front().end(),
                      universialCompareFun<CompareByCongestion>);
            size_t idx=0;
            while(selected.size()<this->_option.populationSize) {
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
                this->_population.erase(i.iterator);
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
