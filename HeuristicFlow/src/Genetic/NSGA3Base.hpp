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

#ifndef Heu_NSGA3BASE_HPP
#define Heu_NSGA3BASE_HPP

#include "NSGA3Abstract.hpp"

namespace Heu {

/**
 * @brief Layers of Reference points
 * 
 */
enum ReferencePointOption {
    SINGLE_LAYER,
    DOUBLE_LAYER
};

template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        RecordOption rOpt,
        PFOption pfOpt,
        ReferencePointOption rpOpt,
        class Args_t,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::mutateFun _mFun_>
class NSGA3Base
    : public NSGA3Abstract<Var_t,ObjNum,DVO,rOpt,pfOpt,Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    NSGA3Base() {};
    virtual ~NSGA3Base() {};
    using Base_t = NSGA3Abstract<Var_t,ObjNum,DVO,rOpt,pfOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_NSGA3ABSTRACT_TYPES

    size_t referencePointPrecision() const {
        return _precision;
    }

    void setReferencePointPrecision(size_t p) {
        _precision=p;
    }

    size_t referencePointCount() const {
        return NchooseK(_precision+this->objectiveNum()-1,_precision);
    }

protected:
    size_t _precision;

    void makeReferencePoses() {
        this->referencePoses.resize(this->objectiveNum(),referencePointCount());
        std::vector<Fitness_t> rfP;
        this->computeReferencePointPoses(this->objectiveNum(),_precision,&rfP);
        std::shuffle(rfP.begin(),rfP.end(),global_mt19937);
        for(size_t c=0;c<this->referencePoses.cols();c++) {
            for(size_t r=0;r<this->referencePoses.rows();r++) {
                this->referencePoses(r,c)=rfP[c][r];
            }
        }
    }
};



template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        RecordOption rOpt,
        PFOption pfOpt,
        class Args_t,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::mutateFun _mFun_>
class NSGA3Base<Var_t,ObjNum,DVO,rOpt,pfOpt,DOUBLE_LAYER,Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
    : public NSGA3Abstract<Var_t,ObjNum,DVO,rOpt,pfOpt,Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    NSGA3Base() {
        _innerPrecision=3;
        _outerPrecision=4;
    };
    virtual ~NSGA3Base() {};
    using Base_t = NSGA3Abstract<Var_t,ObjNum,DVO,rOpt,pfOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_NSGA3ABSTRACT_TYPES

    size_t innerPrecision() const {
        return _innerPrecision;
    }
    
    size_t outerPrecision() const {
        return _outerPrecision;
    }

    void setReferencePointPrecision(size_t i,size_t o) {
        assert(i>=2);
        assert(o>=2);
        _innerPrecision=i;
        _outerPrecision=o;
    }

    size_t referencePointCount() const {
        return NchooseK(_innerPrecision+this->objectiveNum()-1,_innerPrecision)
            +NchooseK(_outerPrecision+this->objectiveNum()-1,_outerPrecision);
    }

protected:
    size_t _innerPrecision;
    size_t _outerPrecision;

    void makeReferencePoses() {
        this->referencePoses.resize(this->objectiveNum(),referencePointCount());
        std::vector<Fitness_t> irfP,orfP;        
        std::shuffle(irfP.begin(),irfP.end(),global_mt19937());
        std::shuffle(orfP.begin(),orfP.end(),global_mt19937());
        this->computeReferencePointPoses(this->objectiveNum(),_innerPrecision,&irfP);
        this->computeReferencePointPoses(this->objectiveNum(),_outerPrecision,&orfP);

        for(int c=0;c<this->referencePoses.cols();c++) {
            for(int r=0;r<this->objectiveNum();r++) {
                if(c<irfP.size()) {
                    this->referencePoses(r,c)=irfP[c][r]*M_SQRT1_2;
                }
                else {
                    this->referencePoses(r,c)=orfP[c-irfP.size()][r];
                }
            }
        }

    }   //  makeReferencePoses()

};

}



#endif  //  Heu_NSGA3BASE_HPP
