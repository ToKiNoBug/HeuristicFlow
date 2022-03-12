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

#ifndef Heu_NSGA2_HPP
#define Heu_NSGA2_HPP

#include "./NSGA2Base.hpp"
#include <type_traits>
namespace Heu
{

/**
   *  @brief NSGA2 MOGA solver. Suitable for not too many objectives.
   *
   *  @tparam Var_t  Type of decisition variable.
   *  @tparam ObjNum Numbers of objectives.
   *  @tparam fOpt Whether greater fitness value means better.
   *  @tparam rOpt Whether the solver records fitness changelog.
   *  @tparam pfOpt Whether to protect the Pareto front from mutation.
   *  @tparam Args_t Type of other parameters.
  */
template<typename Var_t,
         size_t ObjNum,
         DoubleVectorOption DVO=
#ifndef EIGEN_CORE_H
         DoubleVectorOption::Std,
#else
         DoubleVectorOption::Eigen,
#endif
         FitnessOption isGreaterBetter=FITNESS_LESS_BETTER,
         RecordOption Record=DONT_RECORD_FITNESS,
         PFOption ProtectPF=PARETO_FRONT_CAN_MUTATE,
         class Args_t=void,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::initializeFun _iFun_=nullptr,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::fitnessFun _fFun_=nullptr,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::crossoverFun _cFun_=nullptr,
         typename GAAbstract<Var_t,FitnessVec_t<DVO,ObjNum>,Args_t>::mutateFun _mFun_=nullptr>
class NSGA2
    : public NSGA2Base<Var_t,
                    ObjNum,
                    DVO,
                    isGreaterBetter,
                    Record,
                    ProtectPF,
                    Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    using Base_t = NSGA2Base<Var_t,
                    ObjNum,
                    DVO,
                    isGreaterBetter,
                    Record,
                    ProtectPF,
                    Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    using infoUnit2 = typename Base_t::infoUnit2;
    Heu_MAKE_NSGABASE_TYPES

    NSGA2() {
    };

    virtual ~NSGA2() {};

protected:

private:
};

#ifdef EIGEN_CORE_H

/**
 * @brief Partial specialization for NSGA2 using Eigen's array
 * 
 * @tparam Var_t 
 * @tparam ObjNum 
 * @tparam isGreaterBetter 
 * @tparam Record 
 * @tparam ProtectPF 
 * @tparam Args 
 */
template<typename Var_t,
         size_t ObjNum,
         FitnessOption isGreaterBetter,
         RecordOption Record,
         PFOption ProtectPF,
         class Args_t,
         typename GAAbstract<Var_t,FitnessVec_t<DoubleVectorOption::Eigen,ObjNum>,Args_t>::initializeFun _iFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DoubleVectorOption::Eigen,ObjNum>,Args_t>::fitnessFun _fFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DoubleVectorOption::Eigen,ObjNum>,Args_t>::crossoverFun _cFun_,
         typename GAAbstract<Var_t,FitnessVec_t<DoubleVectorOption::Eigen,ObjNum>,Args_t>::mutateFun _mFun_>
class NSGA2<Var_t,
            ObjNum,
            DoubleVectorOption::Eigen,
            isGreaterBetter,
            Record,
            ProtectPF,
            Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
    : public NSGA2Base<Var_t,
                    ObjNum,
                    DoubleVectorOption::Eigen,
                    isGreaterBetter,
                    Record,
                    ProtectPF,
                    Args_t,
            _iFun_,_fFun_,_cFun_,_mFun_>
{
public:
    using Base_t =NSGA2Base<Var_t,
                    ObjNum,
                    DoubleVectorOption::Eigen,
                    isGreaterBetter,
                    Record,
                    ProtectPF,
                    Args_t,
        _iFun_,_fFun_,_cFun_,_mFun_>;
    Heu_MAKE_NSGABASE_TYPES

    using congestComposeFun = typename Base_t::congestComposeFun;
    using infoUnit2 = typename Base_t::infoUnit2;

    NSGA2() {
        setCongestComposeFun();
    };

    virtual ~NSGA2() {};

    inline void setCongestComposeFun(congestComposeFun __ccFun=default_ccFun_liner) {
            this->_ccFun=__ccFun;
    }

    static double default_ccFun_liner(const Fitness_t * f) {
        return f->sum();
    };

    static double default_ccFun_sphere(const Fitness_t * f) {
        return std::sqrt(f->square().sum());
    }

    static double default_ccFun_max(const Fitness_t * f) {
        return f->maxCoeff();
    }

    template<int64_t p>
    static double default_ccFun_powered(const Fitness_t * f) {
        return std::pow(f->power(p).sum(),1.0/p);
    }

    inline void initializePop() {
        if(this->_ccFun==nullptr) {
            setCongestComposeFun();
        }
        Base_t::initializePop();
    }

    virtual Fitness_t bestFitness() const {
        Fitness_t best=Base_t::_population.front()._Fitness;
        for(const Gene & i : Base_t::_population) {
            if(isGreaterBetter) {
                best=best.max(i._Fitness);
            } else {
                best=best.min(i._Fitness);
            }
        }
        return best;
    }

protected:

private:

};

#endif  // Heu_GENETIC_USE_EIGEN
}

#endif // NSGA2_HPP
