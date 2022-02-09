#ifndef OptimT_NSGA3BASE_HPP
#define OptimT_NSGA3BASE_HPP

#include "NSGA3Abstract.hpp"

namespace OptimT {

enum ReferencePointOption {
    SingleLayer,
    DoubleLayer
};

template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        ReferencePointOption rpOpt,
        class Args_t>
class NSGA3Base
    : public NSGA3Abstract<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,Args_t>
{
public:
    NSGA3Abstract() {};
    virtual ~NSGA3Abstract() {};
    using Base_t = NSGA3Abstract<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,Args_t>;
    OptimT_MAKE_NSGA3ABSTRACT_TYPES

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
        this->computeReferencePoses(this->objectiveNum(),_precision,&rfP);
        for(size_t c=0;c<referencePoses.cols();c++) {
            for(size_t r=0;r<referencePoses.rows();r++) {
                this->referencePoses(r,c)=rfP[c][r];
            }
        }
    }
};



template<typename Var_t,
        size_t ObjNum,
        DoubleVectorOption DVO,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class Args_t>
class NSGA3Base<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,DoubleLayer,Args_t>
    : public NSGAAbstract<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,Args_t>
{
public:
    NSGA3Abstract() {
        _innerPrecision=3;
        _outerPrecision=4;
    };
    virtual ~NSGA3Abstract() {};
    using Base_t = NSGAAbstract<Var_t,ObjNum,DVO,fOpt,rOpt,pfOpt,Args_t>;
    OptimT_MAKE_NSGA3ABSTRACT_TYPES

    size_t innerPrecision() const {
        return _innerPrecision;
    }
    
    size_t outerPrecision() const {
        return _outerPrecision
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
        this->computeReferencePoses(this->objectiveNum(),_innerPrecision,&irfP);
        this->computeReferencePoses(this->objectiveNum(),_outerPrecision,&orfP);

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



#endif  //  OptimT_NSGA3BASE_HPP