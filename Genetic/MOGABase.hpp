#ifndef MOGABASE_HPP
#define MOGABASE_HPP

#include "./MOGAAbstract.hpp"

namespace OptimT {

/**
 * @brief MOGA solver base class with compile-time objective count
 * 
 * @tparam Var_t 
 * @tparam ObjNum 
 * @tparam Fitness_t 
 * @tparam fOpt 
 * @tparam rOpt 
 * @tparam pfOpt 
 * @tparam Args 
 */
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

    inline size_t objectiveNum() const {
    return ObjNum;
    }
};


/**
 * @brief MOGA solver base class with runtime objective count
 * 
 * @tparam Var_t 
 * @tparam Fitness_t 
 * @tparam fOpt 
 * @tparam rOpt 
 * @tparam pfOpt 
 * @tparam Args 
 */
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

    inline size_t objectiveNum() const {
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
