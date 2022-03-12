#ifndef Heu_SOGABOX_HPP
#define Heu_SOGABOX_HPP

#include "../EAGlobal/BoxConstraint.hpp"
#include "../EAGlobal/BooleanBox.hpp"
#include "SOGA.hpp"

namespace Heu
{

template<typename VarEle,
    size_t VARDIM,
    DoubleVectorOption DVO,
    BoxShape BS=SQUARE_BOX,
    FitnessOption isGreaterBetter=FITNESS_LESS_BETTER,
    RecordOption Record=DONT_RECORD_FITNESS,
    class otherArgs_t=void>
class SOGABox
    : public SOGA<typename BoxConstraint<VarEle,VARDIM,DVO,otherArgs_t,BS>::Var_t,
        isGreaterBetter,
        Record,
        BoxConstraint<VarEle,VARDIM,DVO,otherArgs_t,BS>>
{
public:
    using Base_t = SOGA<typename BoxConstraint<VarEle,VARDIM,DVO,otherArgs_t,BS>::Var_t,
        isGreaterBetter,
        Record,
        BoxConstraint<VarEle,VARDIM,DVO,otherArgs_t,BS>>;
    
    Heu_MAKE_GABASE_TYPES

    using Box_t = HeuPrivate::pri_BoxConstraint<VarEle,VARDIM,DVO,BS>;

    inline Box_t & box() {
        return this->_args;
    }

    inline const Box_t & box() const {
        return this->_args;
    }

};

}   //  namespace Heu

#endif  //  Heu_SOGABOX_HPP