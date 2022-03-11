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

#ifndef Heu_LAB4NSGA3_H
#define Heu_LAB4NSGA3_H

#include "includes.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>

Eigen::ArrayXd sample2Intercept(Eigen::MatrixXd);
std::vector<Eigen::ArrayXd> makeReferencePoints(const uint64_t dimN,const uint64_t precision);

static const size_t VarDim=4;
static const size_t ObjNum=3;
//static const double alpha=2;
/**
 * @brief DTLZ7
 * 
 */
class testNSGA3
    : public Heu::NSGABase<Eigen::Array<double,VarDim,1>,
        ObjNum,
        Heu::Eigen,
        Heu::FITNESS_LESS_BETTER,
        Heu::RECORD_FITNESS,
        Heu::PARETO_FRONT_CAN_MUTATE,void,nullptr,nullptr,nullptr,nullptr>
{
public:
    testNSGA3() {
        _innerPrecision=2;
        _outerPrecision=3;
    }

    virtual ~testNSGA3() {

    }

    using Base_t = Heu::NSGABase<Eigen::Array<double,VarDim,1>,
    ObjNum,
    Heu::Eigen,
    Heu::FITNESS_LESS_BETTER,
    Heu::RECORD_FITNESS,
    Heu::PARETO_FRONT_CAN_MUTATE,void,nullptr,nullptr,nullptr,nullptr>;

    Heu_MAKE_NSGABASE_TYPES
    using RefPoint = size_t;

    struct infoUnit3 : public infoUnitBase_t
    {
    public:
        Eigen::Array<double,ObjNum,1> translatedFitness;
        size_t closestIdx;
        double distance;
    };

    inline size_t innerPrecision() const {
        return _innerPrecision;
    }

    inline size_t outerPrecisioin() const {
        return _outerPrecision;
    }

    inline void setPrecision(size_t inner,size_t outer) {
        _innerPrecision=inner;
        _outerPrecision=outer;
    }

    Eigen::Array<double,ObjNum,1> bestFitness() const {
        Eigen::Array<double,ObjNum,1> b=this->_population.front()._Fitness;
        for(auto i : this->_population) {
            b=b.min(i._Fitness);
        }
        return b;
    }

    const Eigen::Array<double,ObjNum,Eigen::Dynamic> & referencePoints() const {
        return referencePoses;
    }

public:

    static void iFun(Eigen::Array<double,VarDim,1> * v) {
        v->setRandom();
        (*v)=(*v+1)/2;
    }

    static void fFun(const Eigen::Array<double,VarDim,1> * v,
        Eigen::Array<double,ObjNum,1> * f) {
        f->segment<ObjNum-1>(0)=v->segment<ObjNum-1>(0);
        const double g=1+9/(std::sqrt(f->square().sum())+1e-40)*(f->sum());
        auto f_i=f->segment<ObjNum-1>(0);
        auto a=f_i/(1+g)*(1+(3*M_PI*f_i).sin());
        const double h=ObjNum-a.sum();
        f->operator[](ObjNum-1)=(1+g)*h;        
    }

    static void cFun(const Eigen::Array<double,VarDim,1> * p1,
        const Eigen::Array<double,VarDim,1> * p2,
        Eigen::Array<double,VarDim,1> * c1,
        Eigen::Array<double,VarDim,1> * c2) {
        static const double r=0.2;
        *c1=r*(*p1)+(1-r)*(*p2);
        *c2=r*(*p2)+(1-r)*(*p1);
    }

    static void mFun(Eigen::Array<double,VarDim,1> * v) {
        double & x =v->operator[](size_t(Heu::randD(0,VarDim)));
        x+=Heu::randD(-1,1)*0.05;
        if(x<0)
            x=0;
        if(x>1)
            x=1;
    }

protected:
    size_t _innerPrecision;
    size_t _outerPrecision;
    Eigen::Array<double,ObjNum,Eigen::Dynamic> referencePoses;

    void customOptAfterInitialization() {
        auto inner=makeReferencePoints(ObjNum,_innerPrecision);
        auto outer=makeReferencePoints(ObjNum,_outerPrecision);
        referencePoses.resize(ObjNum,inner.size()+outer.size());
        
        for(size_t c=0;c<referencePoses.cols();c++) {
            if(c<outer.size()) {
                referencePoses.col(c)=outer[c];
            }
            else {
                referencePoses.col(c)=inner[c-outer.size()]*M_SQRT1_2;
            }
        }
    }

    static bool sortByDominatedNum(const infoUnit3 * A,const infoUnit3* B) {
        if(A==B)
            return false;
        return (A->domainedByNum)<(B->domainedByNum);
    }

    void select() {
        const size_t popSizeBef=this->_population.size();
        std::vector<infoUnit3> pop;
        pop.reserve(popSizeBef);
        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            pop.emplace_back();
            pop.back().iterator=it;
            pop.back().closestIdx=-1;
        }

        this->sortSpace.resize(popSizeBef);
        for(size_t i=0;i<popSizeBef;i++) {
            this->sortSpace[i]=pop.data()+i;
        }

        this->calculateDominatedNum();
        this->divideLayers();

        const size_t PFSize=this->pfLayers.front().size();

        if(PFSize<=this->_option.populationSize)
            this->updatePF((const infoUnitBase_t **)this->pfLayers.front().data(),
                           this->pfLayers.front().size());

        std::unordered_set<infoUnit3*> selected;
        selected.reserve(this->_option.populationSize);

        std::vector<infoUnit3*> * FlPtr=nullptr;
        bool needRefPoint=false;

        while(true) {
            if(selected.size()==this->_option.populationSize) {
                needRefPoint=false;
                break;
            }
            if(selected.size()+this->pfLayers.front().size()>this->_option.populationSize) {
                needRefPoint=true;
                FlPtr=(typeof(FlPtr))&this->pfLayers.front();
                break;
            }

            for(infoUnitBase_t* i : this->pfLayers.front()) {
                selected.emplace((infoUnit3*)i);
            }
            this->pfLayers.pop_front();

        }

        if(needRefPoint) {
            ///Normalize procedure
            std::unordered_set<infoUnit3*> Fl;
            Fl.reserve(FlPtr->size());
            for(auto i : *FlPtr) {
                Fl.emplace(i);
            }
            normalize(selected,Fl);
            std::unordered_map<size_t,RefPoint> refPoints;
            refPoints.reserve(referencePoses.cols());
            for(size_t i=0;i<referencePoses.cols();i++) {
                refPoints[i]=0;
            }
            std::unordered_multimap<size_t,infoUnit3*> refPoint2Gene;
            refPoint2Gene.reserve(Fl.size()+selected.size());

            ///Associate procedure
            associate(selected,nullptr);
            associate(Fl,&refPoint2Gene);
            nichePreservation(&selected,&Fl,&refPoints);
        }

        //erase all unselected genes
        for(auto i : this->sortSpace) {
            if(selected.find((infoUnit3*)i)==selected.end()) {
                this->_population.erase(i->iterator);
            }
        }

    }   //  end select

    void normalize(const std::unordered_set<infoUnit3*> & selected,
        const std::unordered_set<infoUnit3*> & Fl) const {
        auto it=(selected.begin());
        Eigen::Array<double,ObjNum,1> ideal=(*it)->iterator->_Fitness;
        Eigen::Array<double,ObjNum,ObjNum> extremePoints;
        std::array<const infoUnit3*,ObjNum> extremePtr;
        extremePtr.fill(*it);
        extremePoints.colwise()=ideal;
        for(auto i : selected) {
            ideal=ideal.min(i->iterator->_Fitness);
            for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
                if(i->iterator->_Fitness[objIdx]>extremePoints(objIdx,objIdx)) {
                    extremePoints.col(objIdx)=i->iterator->_Fitness;
                    extremePtr[objIdx]=i;
                }
            }
        }
        for(auto i : Fl) {
            ideal=ideal.min(i->iterator->_Fitness);
            for(size_t objIdx=0;objIdx<ObjNum;objIdx++) {
                if(i->iterator->_Fitness[objIdx]>extremePoints(objIdx,objIdx)) {
                    extremePoints.col(objIdx)=i->iterator->_Fitness;
                    extremePtr[objIdx]=i;
                }
            }
        }

        bool isSingular;

        {
            std::unordered_set<const void *> set;
            set.reserve(ObjNum);
            for(auto i : extremePtr) {
                set.emplace(i);
            }
            isSingular=(set.size()<ObjNum);
        }

        extremePoints.colwise()-=ideal;

        Eigen::Array<double,ObjNum,1> intercept;
        
        
        if(isSingular) {
            std::cout<<"singular determinat = "<<extremePoints.matrix().determinant()<<std::endl;
            for(size_t r=0;r<ObjNum;r++) {
                intercept[r]=extremePoints(r,r);
            }
        }
        else {
            extremePoints2Intercept(extremePoints,&intercept);
        }


        for(auto i : selected) {
            i->translatedFitness=(i->iterator->_Fitness-ideal)/intercept;
        }
        for(auto i : Fl) {
            i->translatedFitness=(i->iterator->_Fitness-ideal)/intercept;
        }
        
    }

    void associate(const std::unordered_set<infoUnit3*> & st,std::unordered_multimap<size_t,infoUnit3*>* LUT) const {

        for(auto i : st) {
            const auto & w=referencePoses;
            const auto & s=i->translatedFitness;

            auto wT_s=w.matrix().transpose()*s.matrix();
            auto wT_s_w=w.rowwise()*(wT_s.array().transpose());
            Eigen::Array<double,ObjNum,Eigen::Dynamic> norm_wTsw=wT_s_w.rowwise()/(w.colwise().squaredNorm());
            auto s_sub_norm_wTsw=norm_wTsw.colwise()-s;
            auto distance=s_sub_norm_wTsw.colwise().squaredNorm();

            int minDistanceIdx;
            i->distance=distance.minCoeff(&minDistanceIdx);
            i->closestIdx=minDistanceIdx;
            if(LUT!=nullptr)
                LUT->emplace(i->closestIdx,i);
        }
    }

    void nichePreservation(std::unordered_set<infoUnit3*> * selected,
            std::unordered_set<infoUnit3*> * Fl,
            std::unordered_map<size_t,RefPoint> * refPoints) {
        for(auto i : *selected) {
            refPoints->operator[](i->closestIdx)++;
        }
        std::vector<std::unordered_map<size_t,RefPoint>::iterator> minNicheIterators;
        minNicheIterators.reserve(refPoints->size());

        std::vector<infoUnit3*> associatedGenesInFl;
        associatedGenesInFl.reserve(Fl->size());

        while(selected->size()<=this->_option.populationSize) {
            findMinSet(*refPoints,&minNicheIterators);
            auto curRefPoint=minNicheIterators[size_t(Heu::randD(0,minNicheIterators.size()))];
            size_t rhoJ=curRefPoint->second;
            findAssociated(*Fl,curRefPoint,&associatedGenesInFl);


            if(!associatedGenesInFl.empty()) {
                infoUnit3 * pickedGene=nullptr;
                if(rhoJ==0) {
                    //find element in associatedGenesInFl with minimum distance
                    infoUnit3 * minGene=associatedGenesInFl.front();
                    for(auto i : associatedGenesInFl) {
                        if(i->distance<minGene->distance) {
                            minGene=i;
                        }
                    }
                    pickedGene=minGene;
                }
                else {
                    //pick a random member in associatedGenesInFl
                    pickedGene=associatedGenesInFl[size_t(Heu::randD(0,associatedGenesInFl.size()))];
                }
                selected->emplace(pickedGene);
                Fl->erase(pickedGene);
                curRefPoint->second++;
            }
            else {
                refPoints->erase(curRefPoint);
            }
        }   //  end while
    }

    inline static void findMinSet(std::unordered_map<size_t,RefPoint> & refPoints,
            std::vector<std::unordered_map<size_t,RefPoint>::iterator> * minNicheIterators) {
        minNicheIterators->clear();
        size_t minNiche=-1;
        for(auto i : refPoints) {
            minNiche=std::min(minNiche,i.second);
        }
        for(auto it=refPoints.begin();it!=refPoints.end();++it) {
            if(it->second==minNiche) {
                minNicheIterators->emplace_back(it);
            }
        }
    }

    inline static void findAssociated(const std::unordered_set<infoUnit3*> & Fl,
            const std::unordered_map<size_t,RefPoint>::iterator & refP,
            std::vector<infoUnit3*> * associatedGenesInFl) {
        associatedGenesInFl->clear();

        for(auto i : Fl) {
            if(i->closestIdx==refP->first) {
                associatedGenesInFl->emplace_back(i);
            }
        }

    }

    inline static void extremePoints2Intercept(const Eigen::Array<double,ObjNum,ObjNum> & P,
        Eigen::Array<double,ObjNum,1> * intercept) {
        auto P_transpose_inv=P.transpose().matrix().inverse();
        auto ONE=Eigen::Matrix<double,Eigen::Dynamic,1>::Ones(P.cols(),1);
        auto one_div_intercept=(P_transpose_inv*ONE).array();
        *intercept=1.0/one_div_intercept;
    }
};


void testNSGA3Expri();

#endif  //  Heu_LAB4NSGA3_H
