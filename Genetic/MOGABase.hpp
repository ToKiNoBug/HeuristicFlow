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

#ifndef Heu_MOGABASE_HPP
#define Heu_MOGABASE_HPP

#include "./MOGAAbstract.hpp"

#ifndef Heu_NO_RTASSERT
#include <assert.h>
#endif

#ifndef Heu_MOGA_MaxRunTimeObjNum
#define Heu_MOGA_MaxRunTimeObjNum 32
#endif

namespace Heu {

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
        class Args_t>
class MOGABase
        : public MOGAAbstract<Var_t,ObjNum,Fitness_t,fOpt,rOpt,pfOpt,Args_t>
{
public:
    using Base_t = MOGAAbstract<Var_t,ObjNum,Fitness_t,fOpt,rOpt,pfOpt,Args_t>;
    Heu_MAKE_GABASE_TYPES

    MOGABase() {};
    virtual ~MOGABase() {};

    inline constexpr size_t objectiveNum() const {
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
        class Args_t>
class MOGABase<Var_t,Dynamic,Fitness_t,fOpt,rOpt,pfOpt,Args_t>
        : public MOGAAbstract<Var_t,Dynamic,Fitness_t,fOpt,rOpt,pfOpt,Args_t>
{
public:

    MOGABase() {};
    virtual ~MOGABase() {};

    using Base_t = MOGAAbstract<Var_t,Dynamic,Fitness_t,fOpt,rOpt,pfOpt,Args_t>;
    Heu_MAKE_GABASE_TYPES

    inline size_t objectiveNum() const {
    return _objectiveNum;
    }

    inline void setObjectiveNum(size_t _objNum) {
#ifndef Heu_NO_RTASSERT
        assert(_objNum>1);
        assert(_objNum<=Heu_MOGA_MaxRunTimeObjNum);
#endif
        _objectiveNum=_objNum;
    }

protected:
    size_t _objectiveNum;
};


}


#endif // MOGABASE_HPP
