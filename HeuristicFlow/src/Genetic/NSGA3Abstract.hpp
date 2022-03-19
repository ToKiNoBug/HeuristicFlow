// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef Heu_NSGA3ABSTRACT_HPP
#define Heu_NSGA3ABSTRACT_HPP

#include <unordered_map>
#include <unordered_set>
#include "NSGABase.hpp"

namespace Heu
{

namespace internal
{

template<typename Var_t,
        int ObjNum,
        RecordOption rOpt,
        class Args_t,
         typename GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::mutateFun _mFun_>
class NSGA3Abstract
    : public NSGABase<Var_t,ObjNum,FITNESS_LESS_BETTER,rOpt,Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
{
    using Base_t = NSGABase<Var_t,ObjNum,FITNESS_LESS_BETTER,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
public:
    NSGA3Abstract() {};
    virtual ~NSGA3Abstract() {};
    Heu_MAKE_NSGABASE_TYPES
    using RefPointIdx_t = size_t;

    using RefMat_t=Eigen::Array<double,ObjNum,Eigen::Dynamic>;

    inline const RefMat_t & referencePoints() const {
        return referencePoses;
    }

    struct infoUnit3 : public infoUnitBase_t
    {
        Fitness_t translatedFitness;
        size_t closestRefPoint;
        double distance;
    };

protected:
    RefMat_t referencePoses;

    void computeReferencePointPoses(const size_t dimN,
        const size_t precision,
        std::vector<Fitness_t> * dst) const {
        dst->clear();
        dst->reserve(Heu::NchooseK(dimN+precision-1,precision));

        pri_startRP(dimN,precision,dst);
    }

    void select() {
        const size_t popSizeBef=this->_population.size();
        std::vector<infoUnit3> pop;
        pop.reserve(popSizeBef);
        for(auto it=this->_population.begin();it!=this->_population.end();++it) {
            pop.emplace_back();
            pop.back().iterator=it;
            pop.back().closestRefPoint=-1;
        }

        this->sortSpace.resize(popSizeBef);
        for(size_t i=0;i<popSizeBef;i++) {
            this->sortSpace[i]=pop.data()+i;
        }

        this->calculateDominatedNum();
        this->divideLayers();

        const size_t PFSize=this->pfLayers.front().size();

        if(PFSize<=this->_option.populationSize)
            this->updatePF((const infoUnitBase_t **)this->pfLayers.front().data(),this->pfLayers.front().size());

        std::unordered_set<infoUnit3*> selected;
        selected.reserve(this->_option.populationSize);
        std::vector<infoUnit3*> * FlPtr=nullptr;
        bool needRefPoint;

        while(true) {
            if(selected.size()==this->_option.populationSize) {
                needRefPoint=false;
                break;
            }
            if(selected.size()+this->pfLayers.front().size()>this->_option.populationSize) {
                needRefPoint=true;
                FlPtr=(decltype(FlPtr))&this->pfLayers.front();
                break;
            }

            for(infoUnitBase_t* i : this->pfLayers.front()) {
                selected.emplace((infoUnit3*)i);
            }
            this->pfLayers.pop_front();

        }


        if(needRefPoint) {
            ///Normalize procedure
            std::unordered_multimap<RefPointIdx_t,infoUnit3*> Fl;
            Fl.reserve(FlPtr->size());
            normalize(selected,*FlPtr);
            std::unordered_map<RefPointIdx_t,size_t> refPoints;
            refPoints.reserve(referencePoses.cols());
            for(int i=0;i<referencePoses.cols();i++) {
                refPoints.emplace(i,0);
            }

            ///Associate procedure
            associate(selected);
            associate(*FlPtr,&Fl);

            nichePreservation(&selected,&Fl,&refPoints);
        }
        
        //erase all unselected genes
        for(auto i : this->sortSpace) {
            if(selected.find((infoUnit3*)i)==selected.end()) {
                this->_population.erase(i->iterator);
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

    }   //  end selection

    /// to be reloaded
    virtual void normalize(const std::unordered_set<infoUnit3*> & selected,
        const std::vector<infoUnit3*> & Fl) const {

        const size_t M=this->objectiveNum();
        stdContainer<const infoUnit3*,ObjNum> extremePtrs;
        initializeSize<ObjNum>::template resize<decltype(extremePtrs)>(&extremePtrs,M);
        
        Eigen::Array<double,ObjNum,ObjNum> extremePoints;
        
        Fitness_t ideal,intercepts;
        
        extremePoints.resize(M,M);
        ideal.resize(M);
        intercepts.resize(M);
        
        ideal.setConstant(pinfD);

        for(size_t c=0;c<M;c++) {
            extremePtrs[c]=*(Fl.begin());
            for(size_t r=0;r<M;r++) {
                extremePoints(r,c)=extremePtrs[c]->iterator->_Fitness[r];
            }
        }

        for(auto i : selected) {
            ideal=ideal.min(i->iterator->_Fitness);
            for(size_t objIdx=0;objIdx<M;objIdx++) {
                if(i->iterator->_Fitness[objIdx]>extremePtrs[objIdx]->iterator->_Fitness[objIdx]) {
                    extremePtrs[objIdx]=i;
                }
            }
        }

        for(auto i : Fl) {
            ideal=ideal.min(i->iterator->_Fitness);
            for(size_t objIdx=0;objIdx<M;objIdx++) {
                if(i->iterator->_Fitness[objIdx]>extremePtrs[objIdx]->iterator->_Fitness[objIdx]) {
                    extremePtrs[objIdx]=i;
                }
            }
        }

        for(size_t c=0;c<M;c++) {
            extremePoints.col(c)=extremePtrs[c]->iterator->_Fitness-ideal;
        }

        if(isSingular(extremePoints)) {
            for(size_t r=0;r<M;r++) {
                intercepts[r]=extremePoints(r,r);
            }
        }
        else {
            extremePoints2Intercept(extremePoints,&intercepts);
        }

        for(auto i : selected) {
            i->translatedFitness=(i->iterator->_Fitness-ideal)/intercepts;
        }

        for(auto i : Fl) {
            i->translatedFitness=(i->iterator->_Fitness-ideal)/intercepts;
        }
    }

    size_t findNearest(const Fitness_t & s,double * dist) const {
        
        const auto & w=this->referencePoses;
        auto wT_s=w.matrix().transpose()*s.matrix();
        auto wT_s_w=w.rowwise()*(wT_s.array().transpose());
        Eigen::Array<double,ObjNum,Eigen::Dynamic> norm_wTsw
            =wT_s_w.rowwise()/(w.colwise().squaredNorm());
        auto s_sub_norm_wTsw=norm_wTsw.colwise()-s;
        auto distance=s_sub_norm_wTsw.colwise().squaredNorm();

        int minIdx;
        *dist=distance.minCoeff(&minIdx);
        return minIdx;
    }

    /// to be reloaded
    virtual void associate(const std::unordered_set<infoUnit3*> & selected) const {
        for(auto i : selected) {
            i->closestRefPoint=findNearest(i->translatedFitness,&i->distance);
        }
    }

    virtual void associate(const std::vector<infoUnit3*> & Fl_src,
        std::unordered_multimap<RefPointIdx_t,infoUnit3*> * Fl_dst) const {
            for(auto i : Fl_src) {
                RefPointIdx_t idx=findNearest(i->translatedFitness,&i->distance);
                i->closestRefPoint=idx;
                Fl_dst->emplace(idx,i);
            }
    }

    virtual void nichePreservation(std::unordered_set<infoUnit3*> * selected,
        std::unordered_multimap<RefPointIdx_t,infoUnit3*> * Fl,
        std::unordered_map<RefPointIdx_t,size_t> * refPoints) const {
        
        for(auto i : *selected) {
            refPoints->operator[](i->closestRefPoint)++;
        }

        std::vector<std::unordered_map<RefPointIdx_t,size_t>::iterator> minNicheIterators;
        minNicheIterators.reserve(refPoints->size());

        std::pair<typename std::unordered_multimap<RefPointIdx_t,infoUnit3*>::iterator,
            typename std::unordered_multimap<RefPointIdx_t,infoUnit3*>::iterator>
                associatedGenesInFl;
        
        while(selected->size()<this->_option.populationSize) {
            findMinSet(*refPoints,&minNicheIterators);
            auto curRefPoint=minNicheIterators[randIdx(minNicheIterators.size())];
            size_t rhoJ=curRefPoint->second;

            associatedGenesInFl=Fl->equal_range(curRefPoint->first);

            if(associatedGenesInFl.first!=associatedGenesInFl.second) { //  not empty
                typename std::unordered_multimap<RefPointIdx_t,infoUnit3*>::iterator
                        pickedGene;
                if(rhoJ==0) {
                    //find element in associatedGenesInFl with minimum distance
                    typename std::unordered_multimap<RefPointIdx_t,infoUnit3*>::iterator
                        minGene=associatedGenesInFl.first;
                    for(auto it=associatedGenesInFl.first;it!=associatedGenesInFl.second;++it) {
                        if(it->second->distance<minGene->second->distance) {
                            minGene=it;
                        }
                    }
                    pickedGene=minGene;
                }
                else {
                    //pick a random member in associatedGenesInFl
                    size_t N=0;
                    for(auto it=associatedGenesInFl.first;it!=associatedGenesInFl.second;++it) {
                        N++;
                    }
                    for(auto it=associatedGenesInFl.first;it!=associatedGenesInFl.second;++it) {
                        if(randD()*N<=1) {
                            pickedGene=it;
                            break;
                        }
                        N--;
                    }

                }

                selected->emplace(pickedGene->second);
                Fl->erase(pickedGene);
                curRefPoint->second++;
            }
            else {
                refPoints->erase(curRefPoint);
            }
        }   //  end while
    }

    inline static void findMinSet(std::unordered_map<RefPointIdx_t,size_t> & refPoints,
            std::vector<std::unordered_map<RefPointIdx_t,size_t>::iterator> * minNicheIterators) {
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


private:
    void pri_makeRP(const size_t dimN,const size_t precision,
        const size_t curDim,const size_t curP,const size_t accum,
        Fitness_t * rec,
        std::vector<Fitness_t> * dst) const {

        if(curDim+1>=dimN) {
            rec->operator[](dimN-1)=1.0-double(accum)/precision;
            dst->emplace_back(*rec);
            return;
        }

        for(size_t p=0;p+accum<=precision;p++) {
            if(curDim>=0)
                rec->operator[](curDim)=double(p)/precision;
            pri_makeRP(dimN,precision,curDim+1,p,accum+p,rec,dst);
        }
    }

    void pri_startRP(const size_t dimN,
        const size_t precision,
        std::vector<Fitness_t> * dst) const {

        Fitness_t rec;
#if __cplusplus >=201703L
        if constexpr(ObjNum==Runtime) {
            rec.resize(this->objectiveNum());
        }
#else
        if (ObjNum==Eigen::Dynamic) {
            rec.resize(this->objectiveNum());
        }
#endif

        pri_makeRP(dimN,precision,0,0,0,&rec,dst);
    }
    

    inline static bool isSingular(const Eigen::Array<double,ObjNum,ObjNum> & mat) {
        return std::abs(mat.matrix().determinant())<=1e-10;
    }

    inline static void extremePoints2Intercept(const Eigen::Array<double,ObjNum,ObjNum> & P,
        Fitness_t * intercept) {
        auto P_transpose_inv=P.transpose().matrix().inverse();
        auto ONE=Eigen::Matrix<double,ObjNum,1>::Ones(P.cols(),1);
        auto one_div_intercept=(P_transpose_inv*ONE).array();
        *intercept=1.0/one_div_intercept;
    }

};

#define Heu_MAKE_NSGA3ABSTRACT_TYPES \
Heu_MAKE_NSGABASE_TYPES \
using RefMat_t = typename Base_t::RefMat_t; \
using infoUnit3 = typename Base_t::infoUnit3; \
using RefPointIdx_t = typename Base_t::RefPointIdx_t;

}   //  internal

}   //  Heu

#endif  //  Heu_NSGA3ABSTRACT_HPP
