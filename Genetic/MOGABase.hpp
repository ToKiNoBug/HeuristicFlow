#ifndef MOGABASE_HPP
#define MOGABASE_HPP

#include "./MOGAAbstract.hpp"

namespace OptimT {


template<typename Var_t,
        size_t ObjNum,
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class ...Args>
class MOGABase
        : public MOGAAbstract<Var_t,ObjNum,Fitness_t,fOpt,rOpt,pfOpt,Args...>
{
public:
    using Base_t = MOGAAbstract<Var_t,ObjNum,Fitness_t,fOpt,rOpt,pfOpt,Args...>;
    OptimT_MAKE_GABASE_TYPES

    MOGABase() {};
    virtual ~MOGABase() {};

    size_t objectiveNum() const {
    return ObjNum;
    }
};

template<typename Var_t,
        typename Fitness_t,
        FitnessOption fOpt,
        RecordOption rOpt,
        PFOption pfOpt,
        class ...Args>
class MOGABase<Var_t,Dynamic,Fitness_t,fOpt,rOpt,pfOpt,Args...>
        : public MOGAAbstract<Var_t,Dynamic,Fitness_t,fOpt,rOpt,pfOpt,Args...>
{
public:

    MOGABase() {};
    virtual ~MOGABase() {};

    using Base_t = MOGAAbstract<Var_t,Dynamic,Fitness_t,fOpt,rOpt,pfOpt,Args...>;
    OptimT_MAKE_GABASE_TYPES

    size_t objectiveNum() const {
    return _objectiveNum;
    }

    void setObjectiveNum(size_t on) {
        _objectiveNum=on;
    }

protected:
    size_t _objectiveNum;
};


}


#endif // MOGABASE_HPP
