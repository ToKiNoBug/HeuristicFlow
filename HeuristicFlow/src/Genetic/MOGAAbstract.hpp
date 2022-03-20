// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_MOGAABSTRACT_HPP
#define EIGEN_HEU_MOGAABSTRACT_HPP

#include "./GABase.hpp"
#include "../EAGlobal/Pareto.hpp"
#include <queue>
#include <functional>
#include <unordered_set>

namespace Eigen
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
        int ObjNum,
        FitnessOption fOpt,
        RecordOption rOpt,
        class Args_t,
         typename GAAbstract<Var_t,Eigen::Array<double,ObjNum,1>,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,Eigen::Array<double,ObjNum,1>,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,Eigen::Array<double,ObjNum,1>,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,Eigen::Array<double,ObjNum,1>,Args_t>::mutateFun _mFun_>
class MOGAAbstract
    : public GABase<Var_t,Eigen::Array<double,ObjNum,1>,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>
{
private:
    using Base_t = GABase<Var_t,Eigen::Array<double,ObjNum,1>,rOpt,Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    static_assert(ObjNum>0||ObjNum==Eigen::Dynamic,"Invalid template parameter Dim");
    static_assert(ObjNum!=1,"You assigend 1 objective in multi-objective problems");
public:
    MOGAAbstract()  {};
    virtual ~MOGAAbstract() {};
    EIGEN_HEU_MAKE_GABASE_TYPES(Base_t)
    using Fitness_t = Eigen::Array<double,ObjNum,1>;

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
    /**
     * @brief calculate ideal point
     * 
     * @return Fitness_t ideal point
     */
    virtual Fitness_t bestFitness() const {
        Fitness_t best=this->_population.front()._Fitness;
        for(const Gene & i : this->_population) {
            if(fOpt==FitnessOption::FITNESS_GREATER_BETTER) {
                best=best.max(i._Fitness);
            } else {
                best=best.min(i._Fitness);
            }
        }
        return best;
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
    
};  // MOGAAbstract

}   //  internal

}   //  namespace Eigen

#endif //   EIGEN_HEU_MOGAABSTRACT_HPP
