// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef Heu_MOGAABSTRACT_HPP
#define Heu_MOGAABSTRACT_HPP

#include "./GABase.hpp"
#include "../EAGlobal/Pareto.hpp"
#include <queue>
#include <functional>
#include <unordered_set>

namespace Heu
{

namespace internal
{

/**
   *  @brief Base class for multi-objective genetic algorithm solver.
   *
   *  @tparam Var_t  Type of decisition variable.
   *  @tparam ObjNum Numbers of objectives.
   *  @tparam Fitness_t Type of fitness value.
   *  @tparam fOpt Whether greater fitness value means better.
   *  @tparam rOpt Whether the solver records fitness changelog.
   *  @tparam pfOpt Whether to protect the Pareto front from mutation.
   *  @tparam Args_t Type of other parameters.
  */
template<typename Var_t,
        size_t ObjNum,
        FitnessOption fOpt,
        RecordOption rOpt,
        class Args_t,
         typename GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,EigenVecD_t<ObjNum>,Args_t>::mutateFun _mFun_>
class MOGAAbstract
    : public GABase<Var_t,EigenVecD_t<ObjNum>,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    MOGAAbstract()  {};
    virtual ~MOGAAbstract() {};

    using Base_t = GABase<Var_t,EigenVecD_t<ObjNum>,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_GABASE_TYPES
    using Fitness_t = EigenVecD_t<ObjNum>;

    ///get pareto front in vec
    inline void paretoFront(std::vector<Fitness_t> & front) const {
        front.clear();  
        front.reserve(_pfGenes.size());
        for(const Gene* i : _pfGenes) {
            front.emplace_back(i->_Fitness);
        }
        return;
    }

    inline void paretoFront(std::vector<std::pair<const Var_t*,const Fitness_t*>> & front) const {
        front.clear();
        front.reserve(_pfGenes.size());
        for(const Gene* i : _pfGenes) {
            front.emplace_back(std::make_pair(&(i->self),&(i->_Fitness)));
        }
    }

    inline const std::unordered_set<const Gene*> & pfGenes() const {
        return _pfGenes;
    }

    inline void initializePop() {
        this->prevFrontSize=-1;
        this->_pfGenes.clear();
        this->_pfGenes.reserve(this->_option.populationSize*2);
        Base_t::initializePop();
    }

protected:
    size_t prevFrontSize;
    size_t prevPFCheckSum;
    std::unordered_set<const Gene*> _pfGenes;



    virtual size_t makePFCheckSum() const {
        std::vector<const Gene*> pfvec;
        pfvec.reserve(_pfGenes.size());
        for(auto i : _pfGenes) {
            pfvec.emplace_back(i);
        }
        std::sort(pfvec.begin(),pfvec.end());

        
        size_t checkSum=std::hash<const void*>()(pfvec[0]);
        for(size_t i=1;i<pfvec.size();i++) {
            checkSum^=std::hash<const void*>()(pfvec[i]);
        }
        return checkSum;
    }
    
private:
    static_assert(std::integral_constant<bool,(ObjNum!=1)>::value,
    "HeuristicFlow : You assigned single objective in MOGA");


};  // MOGAAbstract

}   //  internal

}   //  Heu

#endif //   MOGAABSTRACT_HPP
