#include "testNsga2.h"

#include <OptimTemplates/SimpleMatrix>
#include <queue>
using namespace OptimT;

testNsga2::testNsga2()
{

}

void testNsga2::select() {

    std::vector<GeneItPlus_t> pop(0);
    pop.reserve(Base_t::_population.size());
    {
        uint32_t emplaced=0;
        for(auto it=Base_t::_population.begin();it!=Base_t::_population.end();++it) {
            //pop.emplace_back(it);
            pop.emplace_back();
            pop.back().iterator=it;
            pop.back().idx=(emplaced++);
            pop.back().paretoRank=-1;
            pop.back().isSelected=false;
            pop.back().congestion={0,0};
        }
    }

    //(A,B)==true means A strong domain B
    MatrixDynamicSize<bool>
            domainMat(Base_t::_population.size(),Base_t::_population.size());

    for(auto & i : domainMat) {
        i=false;
    }

    for(uint32_t r=0;r<Base_t::_population.size();r++) {
        for(uint32_t c=0;c<Base_t::_population.size();c++) {
            if(r==c)
                continue;
            domainMat(r,c)=isBetter(pop[r].iterator->fitness(),pop[c].iterator->fitness());
        }
    }

    std::vector<uint32_t> domainNum(Base_t::_population.size());
    for(uint32_t r=0;r<domainNum.size();r++) {
        domainNum[r]=0;
        for(uint32_t c=0;c<domainNum.size();c++) {
            domainNum[r]+=domainMat(r,c);
        }
    }

    std::list<std::vector<uint32_t>>paretoLayer_Idx;

    uint32_t layeredN=0;
    while(layeredN<Base_t::_population.size()) {
        paretoLayer_Idx.emplace_back();
        paretoLayer_Idx.back().clear();
        paretoLayer_Idx.back().reserve(Base_t::_population.size());

        for(uint32_t idx=0;idx<Base_t::_population.size();idx++) {
            if(pop[idx].paretoRank>=0) {
                //a gene which has been sorted into a layer shouldn't be sorted again
                continue;
            }
            if(domainNum[idx]+1+layeredN>=Base_t::_population.size()) {
                paretoLayer_Idx.back().emplace_back(idx);
                pop[idx].paretoRank=paretoLayer_Idx.size()-1;
            }
        }
        layeredN+=paretoLayer_Idx.back().size();
    }
    if(paretoLayer_Idx.back().size()<=0) {
        paretoLayer_Idx.pop_back();
    }

    std::queue<uint32_t> selectedQueue;
    bool needCalculateCongestion=false;
    for(const auto & curLayer : paretoLayer_Idx) {
        //don't need to calculate congestion
        if(selectedQueue.size()==Base_t::_option.populationSize){
            needCalculateCongestion=false;
            break;
        }
        //need to calculate congestion
        if(selectedQueue.size()+curLayer.size()>Base_t::_option.populationSize) {
            needCalculateCongestion=true;
            break;
        }
        for(auto i : curLayer) {
            selectedQueue.emplace(i);
        }
    }

    //calculate congestion
    if(needCalculateCongestion) {
        //using Var_t = std::array<double,2>;
        using sortUnit = std::pair<GeneItPlus_t*,uint8_t>;//Var_t pointer and object idx
        std::vector<sortUnit> sortSpace;
        sortSpace.resize(pop.size());
        //make sort space
        for(uint32_t idx=0;idx<pop.size();idx++) {
            sortSpace[idx].first=&(pop[idx]);
        }
        bool (*cmpFun)(const sortUnit & A,const sortUnit & B)
                =[](const sortUnit & A,const sortUnit & B) {
            //compare on single object
            return (A.first->iterator->_Fitness[A.second])
                    <B.first->iterator->_Fitness[A.second];};

        //sort on n-th object
        for(uint8_t objIdx=0;objIdx<2;objIdx++) {
            for(uint32_t idx=0;idx<pop.size();idx++) {
                sortSpace[idx].second=objIdx;
            }
            std::sort(sortSpace.begin(),sortSpace.end(),cmpFun);
            sortSpace.front().first->congestion[objIdx]=1.0/0.0;
            sortSpace.back().first->congestion[objIdx]=1.0/0.0;

            //fill in congestion on single dim
            for(uint32_t idx=1;idx+1<sortSpace.size();idx++) {
                sortSpace[idx].first->congestion[objIdx]=
                        std::abs(
                            sortSpace[idx-1].first->iterator->_Fitness[objIdx]-
                        sortSpace[idx+1].first->iterator->_Fitness[objIdx]);
            }

        }

        //sort through congestion
        using congestSortUnit = std::pair<GeneItPlus_t*,double>;
        std::vector<congestSortUnit> congestSortSpace;
        congestSortSpace.resize(pop.size());

        for(auto & i : congestSortSpace) {
            i.second=0;
            for(double j : i.first->congestion) {
                i.second+=j;
            }
        }

        std::sort(congestSortSpace.begin(),congestSortSpace.end(),
                [](const congestSortUnit & A,const congestSortUnit & B) {
            return A.second<B.second;});

        //fill selected queue
        for(const auto & i : congestSortSpace) {
            if(selectedQueue.size()>=Base_t::_option.populationSize)
                break;
            selectedQueue.emplace(i.first->idx);
        }
    }//congestion finished

    //mark selected as selected
    while(!selectedQueue.empty()) {
        pop[selectedQueue.front()].isSelected=true;
        selectedQueue.pop();
    }

    //erase unselected
    for(const auto & it : pop) {
        if(!it.isSelected) {
            Base_t::_population.erase(it.iterator);
        }
    }
}

bool testNsga2::isBetter(const std::array<double,2>& A,const std::array<double,2>& B) {
    return (A[0]<B[0])&&(A[1]<B[1]);
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

    for(uint32_t r=0;r<Base_t::_population.size();r++) {
        domainedNum[r]=0;
        for(uint32_t c=0;c<Base_t::_population.size();c++) {
            if(r==c) {
                continue;
            }
            domainedNum[r]+=!isBetter(pop[r]->_Fitness,pop[c]->_Fitness);
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

void runNSGA2() {
    testNsga2 algo;
    using Var_t = std::array<double,2>;
    using Args_t = testNsga2::ArgsType;

    //object1: f1=4x^2+4y^2
    //object2: f2=(x-5)^2+(y-5)^2
    //0<=x<=5,0<=y<=3
    algo.initialize(
                [](Var_t* v,const Args_t *) {
        v->at(0)=OtGlobal::randD(0,5);
        v->at(1)=OtGlobal::randD(0,3);},
    [](const Var_t * v,const Args_t *) {
        std::array<double,2> fitness;
        const double & x=v->at(0),y=v->at(1);
        fitness[0]=4*x*x+4*y*y;
        fitness[1]=OT_square(x-5)+OT_square(y-5);
        return fitness;},
    [](const Var_t * p1,const Var_t * p2,Var_t * c1,Var_t * c2,const Args_t *) {
        for(uint32_t idx=0;idx<p1->size();idx++) {
            c1->at(idx)=(OtGlobal::randD()<0.5)?p1->at(idx):p2->at(idx);
            c2->at(idx)=(OtGlobal::randD()<0.5)?p1->at(idx):p2->at(idx);
        } },
    [](Var_t * v,const Args_t *) {
        v->at(0)+=OtGlobal::randD(-0.05,0.05);
        v->at(0)=std::max(v->at(0),0.0);
        v->at(0)=std::min(5.0,v->at(0));
        v->at(1)+=OtGlobal::randD(-0.05,0.05);
        v->at(1)=std::max(v->at(1),0.0);
        v->at(1)=std::min(3.0,v->at(1));
    },
    nullptr,
    GAOption()
    );

    algo.run();

}
