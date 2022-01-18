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

#ifndef NSGA2_HPP
#define NSGA2_HPP

#include "./GABase.hpp"
#include <queue>
#include <unordered_set>
namespace OptimT
{
///whether to protect pareto front when mutation or not
enum PFOption : unsigned char {
    PARETO_FRONT_DONT_MUTATE=true,
    PARETO_FRONT_CAN_MUTATE=false
};

///NSGA2 MOGA solver
template<typename Var_t,size_t ObjNum,
         FitnessOption isGreaterBetter,
         RecordOption Record,
         PFOption ProtectPF,
         class ...Args>
class NSGA2 : public GABase<Var_t,std::array<double,ObjNum>,Record,Args...>
{
public:
    using This_t = NSGA2;
    using Base_t = GABase<Var_t,std::array<double,ObjNum>,Record,Args...>;
    using Fitness_t = std::array<double,ObjNum>;
    OPTIMT_MAKE_GABASE_TYPES

    using congestComposeFun = double(*)(const Fitness_t *,const ArgsType*);

    NSGA2() {
        Base_t::_initializeFun=
                [](Var_t*,const ArgsType*){};
        Base_t::_fitnessFun=
                [](const Var_t*,const ArgsType*,Fitness_t*){};
        Base_t::_crossoverFun=
                [](const Var_t*,const Var_t*,Var_t*,Var_t*,const ArgsType*){};
        Base_t::_mutateFun=
                [](Var_t*,const ArgsType*){};
        Base_t::_otherOptFun=
                [](ArgsType*,std::list<Gene>*,size_t,size_t,const GAOption*){};
        _ccFun=default_ccFun_liner;
    };

    virtual ~NSGA2() {};

    void setCongestComposeFun(congestComposeFun __ccFun=default_ccFun_liner) {
            _ccFun=__ccFun;
    }

    static double default_ccFun_liner(const Fitness_t * f,const ArgsType*) {
        double result=f->at(0);
        for(size_t objIdx=1;objIdx<ObjNum;objIdx++) {
            result+=f->at(objIdx);
        }
        return result;
    };

    static double default_ccFun_sphere(const Fitness_t * f,const ArgsType*) {
        double result=OT_square(f->at(0));
        for(size_t objIdx=1;objIdx<ObjNum;objIdx++) {
            result+=OT_square(f->at(objIdx));
        }
        return std::sqrt(result);
    }

    static double default_ccFun_max(const Fitness_t * f,const ArgsType*) {
        double result=f->at(0);
        for(size_t objIdx=1;objIdx<ObjNum;objIdx++) {
            result=std::max(f->at(objIdx),result);
        }
        return result;
    }

    template<int power>
    static double default_ccFun_powered(const Fitness_t * f,const ArgsType*) {
        double result=0;
        for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
            result=std::pow(f->at(objIdx),(power));
        }
        return std::pow(result,1.0/power);
    }

    virtual Fitness_t bestFitness() const {
        Fitness_t best=Base_t::_population.front()._Fitness;
        for(const Gene & i : Base_t::_population) {
            if(isGreaterBetter) {
                for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
                    best[objIdx]=std::max(best[objIdx],i._Fitness[objIdx]);
                }
            } else {
                for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
                    best[objIdx]=std::min(best[objIdx],i._Fitness[objIdx]);
                }
            }
        }
        return best;
    }

    ///get pareto front in vec
    void paretoFront(std::vector<Fitness_t> & front) const {
        front.clear();  front.reserve(_pfGenes.size());
        for(const Gene* i : _pfGenes) {
            front.emplace_back(i->_Fitness);
        }
        return;
    }

    void paretoFront(std::vector<std::pair<const Var_t*,const Fitness_t*>> & front) const {
        front.clear();
        front.reserve(_pfGenes.size());
        for(const Gene* i : _pfGenes) {
            front.emplace_back(std::make_pair(&(i->self),&(i->_Fitness)));
        }
    }

    const std::unordered_set<const Gene*> & pfGenes() const {
        return _pfGenes;
    }

    ///temporary struct to store infos when selection
    struct infoUnit
    {
    public:
        bool isSelected;
        size_t sortIdx;
        //uint32_t index;
        ///Genes in population that strong domain this gene
        size_t domainedByNum;
        GeneIt_t iterator;
        Fitness_t congestion;
    };

protected:
    size_t prevFrontSize;
    std::unordered_set<const Gene*> _pfGenes;
    congestComposeFun _ccFun;

    virtual void customOptWhenInitialization() {
        prevFrontSize=-1;
        _pfGenes.clear();
        _pfGenes.reserve(Base_t::_option.populationSize*2);
        if(_ccFun==nullptr) {
            setCongestComposeFun();
        }
    }

    ///whether A strong domainates B
    static bool isStrongDomain(const Fitness_t * A,const Fitness_t * B) {
        //if(A==B) return false;
        for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
            if(isGreaterBetter) {
                //if any single fitness of A isn't better than B, A doesn't strong domain B
                if((A->at(objIdx))<(B->at(objIdx))) {
                    return false;
                }
            } else {
                //if any single fitness of A isn't better than B, A doesn't strong domain B
                if((A->at(objIdx))>(B->at(objIdx))) {
                    return false;
                }
            }
        }
        return true;
    } //isStrongDomain

    ///compare by  paretoLayers
    static bool compareFun_DomainedBy(const infoUnit * A,const infoUnit * B) {
        if(A==B) return false;
        return A->domainedByNum<B->domainedByNum;
    }

    ///compare by fitness on single object
    static bool compareFun_Fitness(const infoUnit * A,const infoUnit * B) {
        if(A==B) return false;
        //if(!isGreaterBetter)
            return A->iterator->_Fitness[A->sortIdx]<B->iterator->_Fitness[B->sortIdx];
        //else
            //return A->iterator->_Fitness[A->sortIdx]>B->iterator->_Fitness[B->sortIdx];
    }

    ///compare by congestion
    static bool compareFun_Congestion(const infoUnit * A,const infoUnit * B) {
        return A->congestion[0]>B->congestion[0];
    }

    ///fast nondominated sorting
    virtual void select() {
        const size_t popSizeBefore=Base_t::_population.size();
        std::vector<infoUnit> pop;
        pop.clear();pop.reserve(popSizeBefore);

        for(auto it=Base_t::_population.begin();it!=Base_t::_population.end();++it) {
            pop.emplace_back();
            pop.back().isSelected=false;
            pop.back().sortIdx=0;
            //pop.back().index=pop.size()-1;
            //pop.back().domainedByNum=0;
            pop.back().iterator=it;
        }

        std::vector<infoUnit*> sortSpace(popSizeBefore);
        //make sortspace
        for(size_t i=0;i<popSizeBefore;i++) {
            sortSpace[i]=pop.data()+i;
        }

        //calculate domainedByNum
        for(size_t ed=0;ed<popSizeBefore;ed++) {
            pop[ed].domainedByNum=0;
            for(size_t er=0;er<popSizeBefore;er++) {
                if(er==ed)
                    continue;
                pop[ed].domainedByNum+=
                        isStrongDomain(&(pop[er].iterator->_Fitness),
                                       &(pop[ed].iterator->_Fitness));
            }
        }
        //sort by paretoLayers
        std::sort(sortSpace.begin(),sortSpace.end(),
                compareFun_DomainedBy);

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

        {
            const size_t curFrontSize=paretoLayers.front().size();
            if(prevFrontSize!=curFrontSize) {
                Base_t::_failTimes=0;
                prevFrontSize=curFrontSize;
            } else {
                Base_t::_failTimes++;
            }
            _pfGenes.clear();
            for(const auto i :paretoLayers.front()) {
                _pfGenes.emplace(&*(i->iterator));
            }
        }


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
            for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
                for(infoUnit & i : pop) {
                    i.sortIdx=objIdx;
                }
                std::sort(sortSpace.begin(),sortSpace.end(),compareFun_Fitness);

                sortSpace.front()->congestion[objIdx]=OptimT::pinfD;
                sortSpace.back()->congestion[objIdx]=OptimT::pinfD;

                //calculate congestion on single object


                for(size_t idx=1;idx<popSizeBefore-1;idx++) {

                    sortSpace[idx]->congestion[objIdx]=std::abs(
                                sortSpace[idx-1]->iterator->_Fitness[objIdx]
                               -sortSpace[idx+1]->iterator->_Fitness[objIdx]

                                );
                }
            } // end sort on objIdx

            for(infoUnit & i : pop) {
                //store final congestion at the first congestion
                i.congestion[0]=_ccFun(&i.congestion,&Base_t::args());
                /*
                for(size_t objIdx=1;objIdx<ObjNum;objIdx++) {
                    i.congestion[0]+=i.congestion[objIdx];
                }*/
            }

            std::sort(paretoLayers.front().begin(),paretoLayers.front().end(),
                      compareFun_Congestion);
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

    ///mutate operation
    virtual void mutate() {
        for(auto it=Base_t::_population.begin();it!=Base_t::_population.end();++it) {
            if(OtGlobal::randD()<=Base_t::_option.mutateProb) {
                if(ProtectPF){
                    if(_pfGenes.find(&*it)!=_pfGenes.end()) {
                        continue;
                    }
                }
                Base_t::_mutateFun(&it->self,&Base_t::args());
                it->setUncalculated();
            }
        }
    }


#ifndef OptimT_NO_STATICASSERT
    static_assert(std::integral_constant<bool,(ObjNum>1)>::value,
    "OptimTemplates : You used less than 2 object functions in NSGA2");
#endif
};

}

#endif // NSGA2_HPP
