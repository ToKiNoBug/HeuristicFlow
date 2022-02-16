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

#include "lab4NSGA2.h"
#include <queue>
#include <ctime>
using namespace Heu;

testNsga2::testNsga2() {

}

void testNsga2::select() {

    //if sort by object, this number will be 0 or positive
    static const int32_t
            SORT_BY_DOMAINEDBYNUM=-1,
            SORT_BY_CONGESTION=-2,
            NONE=-2147483648;
    ///popSize before selection
    const size_t popSizeBef=_population.size();
    std::vector<infoUnit> pop;
    pop.clear();
    pop.reserve(popSizeBef);
    for(auto it=_population.begin();it!=_population.end();++it) {
        pop.emplace_back();
        pop.back().index=pop.size()-1;
        //pop.back().isLayered=false;
        pop.back().isSelected=false;
        pop.back().iterator=it;
        pop.back().sortType=NONE;
    }

    //calculate domainedByNum
    for(uint32_t ed=0;ed<popSizeBef;ed++) {
        pop[ed].domainedByNum=0;
        for(uint32_t er=0;er<popSizeBef;er++) {
            if(er==ed) {
                continue;
            }
            pop[ed].domainedByNum
                    +=isBetter(pop[er].iterator->_Fitness,pop[ed].iterator->_Fitness);
        }
    }

    std::vector<const infoUnit *> sortSpace(popSizeBef);
    for(uint32_t idx=0;idx<popSizeBef;idx++) {
        sortSpace[idx]=pop.data()+idx;
        pop[idx].sortType=SORT_BY_DOMAINEDBYNUM;
    }

    bool (*universalCompareFun)(const infoUnit *,const infoUnit *)
            =[](const infoUnit * A,const infoUnit * B) {
        if(A==B) {
            return false;
        }
        //sort on single object
       if(A->sortType>=0) {
           return A->iterator->_Fitness[A->sortType]<B->iterator->_Fitness[B->sortType];
       }
       if(A->sortType==SORT_BY_DOMAINEDBYNUM) {
           return A->domainedByNum<B->domainedByNum;
       }
       //use the first congestion value as total congestion value
       if(A->sortType==SORT_BY_CONGESTION) {
           return A->congestion[0]>B->congestion[0];
       }
       std::cerr<<"SORTING FAULT : UNKNOWN SORTTYPE "<<A->sortType<<std::endl;
       return false;
    };

    //sort by domainedByNum
    std::sort(sortSpace.begin(),sortSpace.end(),universalCompareFun);

    std::list<std::vector<const infoUnit *>> paretoLayers;
    //make by layers
    {
        uint32_t prevDomainedByNum=-1;
        for(const auto i : sortSpace) {
            if(i->domainedByNum!=prevDomainedByNum) {
                paretoLayers.emplace_back();
                paretoLayers.back().clear();
                paretoLayers.back().reserve(popSizeBef);
                prevDomainedByNum=i->domainedByNum;
            }
            paretoLayers.back().emplace_back(i);
        }
    }

    std::queue<const infoUnit *> selected;
    bool needCongestion=true;
    while(true) {
        if(selected.size()==_option.populationSize) {
            needCongestion=false;
            break;
        }
        //must select some from a single layer
        if(selected.size()+paretoLayers.front().size()>_option.populationSize) {
            needCongestion=true;
            break;
        }

        //select them from layers to layers
        for(const auto i : paretoLayers.front()) {
            selected.emplace(i);
        }
        paretoLayers.pop_front();
    }

    //finished selection from a single layer
    if(needCongestion) {
        for(uint32_t obj=0;obj<2;obj++) {
            for(auto & i : pop) {
                i.sortType=obj;
            }

            std::sort(sortSpace.begin(),sortSpace.end(),universalCompareFun);

            pop[sortSpace.front()->index].congestion[obj]=pinfD;
            pop[sortSpace.back()->index].congestion[obj]=pinfD;

            for(uint32_t idx=1;idx<popSizeBef-1;idx++) {
                pop[sortSpace[idx]->index].congestion[obj]=std::abs(
                            pop[sortSpace[idx-1]->index].iterator->_Fitness[obj]
                        -   pop[sortSpace[idx+1]->index].iterator->_Fitness[obj]
                            );
            }
        }

        //summary congestion
        for(auto & i : pop) {
            i.congestion[0]*=1;
            for(uint32_t obj=1;obj<2;obj++) {
                i.congestion[0]+=i.congestion[obj];
            }
            i.sortType=SORT_BY_CONGESTION;
        }

        std::sort(paretoLayers.front().begin(),paretoLayers.front().end(),universalCompareFun);

        for(auto i : paretoLayers.front()) {
            if(selected.size()>=_option.populationSize) {
                break;
            }
            selected.emplace(i);
        }

    }


    //erase unselected
    while(!selected.empty()) {
        pop[selected.front()->index].isSelected=true;
        selected.pop();
    }

    for(auto & i : pop) {
        if(!i.isSelected) {
            _population.erase(i.iterator);
        }
    }

}

bool testNsga2::isBetter(const std::array<double,2>& A,const std::array<double,2>& B) {
    return (A[0]<B[0])&&(A[1]<B[1]);
}

void testNsga2::mutate() {

    for(auto & i : _population) {
        if(Heu::randD()<_option.mutateProb) {
            _mFun(&i.self);
            i.setUncalculated();
        }
    }

}

void testNsga2::paretoFront(std::vector<const Base_t::Gene *> & dst) const {
    std::vector<std::list<Base_t::Gene>::const_iterator> pop;
    pop.clear();
    pop.reserve(Base_t::_population.size());

    for(auto it=Base_t::_population.cbegin();it!=Base_t::_population.cend();++it) {
        pop.emplace_back(it);
    }

    std::vector<uint32_t> domainedNum;
    domainedNum.resize(Base_t::_population.size());

    for(uint32_t ed=0;ed<Base_t::_population.size();ed++) {
        domainedNum[ed]=0;
        for(uint32_t er=0;er<Base_t::_population.size();er++) {
            if(ed==er) {
                continue;
            }
            domainedNum[ed]+=isBetter(pop[er]->_Fitness,pop[ed]->_Fitness);
        }
    }

    std::queue<uint32_t> frontIdx;

    for(uint32_t i=0;i<domainedNum.size();i++) {
        if(domainedNum[i]<=0) {
            frontIdx.emplace(i);
        }
    }

    dst.clear();
    dst.reserve(frontIdx.size());

    while(!frontIdx.empty()) {
        dst.emplace_back(&*pop[frontIdx.front()]);
        frontIdx.pop();
    }

}

std::array<double,2> testNsga2::bestFitness() const {
    std::array<double,2> best=_population.front()._Fitness;
    for(const auto & i : _population) {
        best[0]=std::min(i._Fitness[0],best[0]);
        best[1]=std::min(i._Fitness[1],best[1]);
    }
    return best;
}

void runNSGA2() {
    testNsga2 algo;
    using Var_t = std::array<double,2>;
    using Args_t = testNsga2::ArgsType;

    //object1: f1=4x^2+4y^2
    //object2: f2=(x-5)^2+(y-5)^2
    //0<=x<=5,0<=y<=3

    GAOption opt;
    opt.maxGenerations=3000;
    opt.maxFailTimes=-1;
    opt.populationSize=200;

    algo.setiFun([](Var_t* v) {
        v->at(0)=randD(0,5);
        v->at(1)=randD(0,3);});

    algo.setfFun(    [](const Var_t * v,std::array<double,2>*fitness) {
        const double & x=v->at(0),y=v->at(1);
        fitness->at(0)=4*x*x+4*y*y;
        fitness->at(1)=square(x-5)+square(y-5);
        });

    algo.setcFun(    [](const Var_t * p1,const Var_t * p2,Var_t * c1,Var_t * c2) {
        for(uint32_t idx=0;idx<p1->size();idx++) {
            c1->at(idx)=(randD()<0.5)?p1->at(idx):p2->at(idx);
            c2->at(idx)=(randD()<0.5)?p1->at(idx):p2->at(idx);
        } });

    algo.setmFun(    [](Var_t * v) {
        v->at(0)+=randD(-0.05,0.05);
        v->at(0)=std::max(v->at(0),0.0);
        v->at(0)=std::min(5.0,v->at(0));
        v->at(1)+=randD(-0.05,0.05);
        v->at(1)=std::max(v->at(1),0.0);
        v->at(1)=std::min(3.0,v->at(1));
    });

    algo.setOption(opt);
    algo.initializePop();

    std::cout<<"Start running..."<<std::endl;
    std::clock_t t=std::clock();
    algo.run();
    std::vector<const testNsga2::Gene*> paretoFront;
    algo.paretoFront(paretoFront);
    t=std::clock()-t;
    std::cout<<"Time elapsed:"<<double(t)/CLOCKS_PER_SEC<<std::endl;
    /*
    std::sort(paretoFront.begin(),paretoFront.end(),
              [](const testNsga2::Gene* a,const testNsga2::Gene* b)
    {
        return a->_Fitness[0]<b->_Fitness[0];
    }
              );
    */
    std::cout<<"\n\n\nparetoFront=[";
    for(auto i : paretoFront) {
        std::cout<<i->fitness()[0]<<" , "<<i->fitness()[1]<<";\n";
    }
    std::cout<<"];"<<std::endl;

    std::cout<<"\n\n\n\nHistory=[";
    for(const auto & i : algo.record()) {
        std::cout<<i[0]<<" , "<<i[1]<<";\n";
    }
    std::cout<<"];"<<std::endl;

}
