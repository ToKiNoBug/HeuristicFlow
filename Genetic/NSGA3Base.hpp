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

enum ReferencePointOption {
    SingleLayer,
    DoubleLayer
};

template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        RecordOption rOpt,
        PFOption pfOpt,
        ReferencePointOption rpOpt,
        class Args_t>
class NSGA3Base
    : public NSGA3Abstract<Var_t,ObjNum,DVO,rOpt,pfOpt,Args_t>
{
public:
    NSGA3Base() {};
    virtual ~NSGA3Base() {};
    using Base_t = NSGA3Abstract<Var_t,ObjNum,DVO,rOpt,pfOpt,Args_t>;
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
        class Args_t>
class NSGA3Base<Var_t,ObjNum,DVO,rOpt,pfOpt,DoubleLayer,Args_t>
    : public NSGA3Abstract<Var_t,ObjNum,DVO,rOpt,pfOpt,Args_t>
{
public:
    NSGA3Base() {
        _innerPrecision=3;
        _outerPrecision=4;
    };
    virtual ~NSGA3Base() {};
    using Base_t = NSGA3Abstract<Var_t,ObjNum,DVO,rOpt,pfOpt,Args_t>;
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
            +NchooseK(_outerPrecision,this->objectiveNum()-1,_outerPrecision);
    }

protected:
    size_t _innerPrecision;
    size_t _outerPrecision;

    void makeReferencePoses() {
        this->referencePoses.resize(this->objectiveNum(),referencePointCount());
        std::vector<Fitness_t> irfP,orfP;        
        this->computeReferencePointPoses(this->objectiveNum(),_innerPrecision,&irfP);
        this->computeReferencePointPoses(this->objectiveNum(),_outerPrecision,&orfP);

        for(size_t c=0;c<this->referencePoses.cols();c++) {
            for(size_t r=0;r<this->objectiveNum();r++) {
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
