#ifndef NSGA2_HPP
#define NSGA2_HPP

#include "./GABase.hpp"
#include <queue>
namespace OptimT
{

template<typename Var_t,size_t ObjNum,bool isGreaterBetter,bool Record,class ...Args>
class NSGA2 : public GABase<Var_t,std::array<double,ObjNum>,Record,Args...>
{
public:
    using This_t = NSGA2;
    using Base_t = GABase<Var_t,std::array<double,ObjNum>,Record,Args...>;
    using Fitness_t = std::array<double,ObjNum>;
    OPTIMT_MAKE_GABASE_TYPES

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
    };

    ~NSGA2() {};

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
        const size_t popSize=Base_t::_population.size();
        std::vector<const Gene*> pop;
        pop.clear();    pop.reserve(popSize);
        for(const auto & i : Base_t::_population) {
            pop.emplace_back(&i);
        }
        std::vector<size_t> domainedByNum(popSize);

        for(size_t ed=0;ed<popSize;ed++) {
            domainedByNum[ed]=0;
            for(size_t er=0;er<popSize;er++) {
                if(ed==er)
                    continue;
                domainedByNum[ed]+=
                        isStrongDomain(&pop[er]->_Fitness,&pop[ed]->_Fitness);
            }
        }

        front.clear();  front.reserve(popSize);
        for(size_t i=0;i<popSize;i++) {
            if(domainedByNum[i]<=0) {
                front.emplace_back(pop[i]->_Fitness);
            }
        }
        return;
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
    ///whether A strong domains B
    static bool isStrongDomain(const Fitness_t * A,const Fitness_t * B) {
        for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
            if(isGreaterBetter) {
                //if any single object fitness of A is worse than B, A doesn't strong domain B
                if(A->at(objIdx)<=B->at(objIdx)) {
                    return false;
                }
            } else {
                //if any single object fitness of A is worse than B, A doesn't strong domain B
                if(A->at(objIdx)>=B->at(objIdx)) {
                    return false;
                }

            }
        }
        return true;
    } //isStrongDomain

    ///compare by  paretoLayers
    static bool compareFun_DomainedBy(const infoUnit * A,const infoUnit * B) {
        return A->domainedByNum<B->domainedByNum;
    }

    ///compare by fitness on single object
    static bool compareFun_Fitness(const infoUnit * A,const infoUnit * B) {
        //if(!isGreaterBetter)
            return A->iterator->_Fitness[A->sortIdx]<B->iterator->_Fitness[B->sortIdx];
        //else
            //return A->iterator->_Fitness[A->sortIdx]>B->iterator->_Fitness[B->sortIdx];
    }

    ///compare by congestion
    static bool compareFun_Congestion(const infoUnit * A,const infoUnit * B) {
        return A->congestion[0]>B->congestion[0];
    }

    ///fast non-domain sorting
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
                for(size_t objIdx=1;objIdx<ObjNum;objIdx++) {
                    i.congestion[0]+=i.congestion[objIdx];
                }
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
                Base_t::_mutateFun(&it->self,&Base_t::args());
                it->setUncalculated();
            }
        }
    }

};

}

#endif // NSGA2_HPP
