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
    NSGA3Abstract() {};
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

};

}



#endif  //  OptimT_NSGA3BASE_HPP