// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_PSO_HPP
#define Heu_PSO_HPP
#include "PSOOption.hpp"
#include "PSOBase.hpp"
#include <array>
#include <vector>
#include <tuple>
#include <type_traits>

namespace Heu {

///Generalized PSO solver
/**
 * @brief Generalized PSO solver.
 * 
 * @tparam Var_t Type of determination vector.
 * @tparam DIM Dimensional of Var_t
 * @tparam FitnessOpt Trainning direction
 * @tparam RecordOpt Record trainning curve or not.
 * @tparam Arg_t Any other parameters.
 */
template<class Var_t,
         size_t DIM,
         FitnessOption FitnessOpt=FITNESS_LESS_BETTER,
         RecordOption RecordOpt=DONT_RECORD_FITNESS,
         class Arg_t=void,
         typename internal::PSOParameterPack<Var_t,double,Arg_t>::iFun_t _iFun_=nullptr,
         typename internal::PSOParameterPack<Var_t,double,Arg_t>::fFun_t _fFun_=nullptr>
class PSO : public internal::PSOBase<Var_t,DIM,double,RecordOpt,Arg_t,_iFun_,_fFun_>
{
public:
    using Base_t = internal::PSOBase<Var_t,DIM,double,RecordOpt,Arg_t,_iFun_,_fFun_>;
    Heu_MAKE_PSOABSTRACT_TYPES

    static const DoubleVectorOption Flag =
        (std::is_same<Var_t,stdVecD_t<DIM>>::value)?
        DoubleVectorOption::Std : DoubleVectorOption::Custom;
    
    void setPVRange(double pMin,double pMax,double vMax) {
        for(size_t i=0;i<this->dimensions();i++) {
            this->_posMin[i]=pMin;
            this->_posMax[i]=pMax;
            this->_velocityMax[i]=vMax;
        }
    }
    
    ///function used to provide a result for recording
    virtual double bestFitness() const {
        return this->gBest.fitness;
    }

    ///Candidate function for initialization
    static void default_iFun(Var_t *x,Var_t *v,
    const Var_t * xMin,const Var_t * xMax,
    const Var_t *,const Args_t*) {
        for(size_t idx=0;idx<xMin->size();idx++) {
            x->operator[](idx)=randD(xMin->operator[](idx),xMax->operator[](idx));
            v->operator[](idx)=0;
        }
    }

protected:
    static bool isBetterThan(double a,double b) {
        if(FitnessOpt==FitnessOption::FITNESS_GREATER_BETTER) {
            return a>b;
        } else {
            return a<b;
        }
    }

    virtual void updatePGBest() {
        Point_t * curGBest=&this->_population.front().pBest;

        for(Particle_t & i : this->_population) {
            if(isBetterThan(i.fitness,i.pBest.fitness)) {
                i.pBest=i;
            }

            if(isBetterThan(i.pBest.fitness,curGBest->fitness)) {
                curGBest=&i.pBest;
            }
        }

            if(isBetterThan(curGBest->fitness,this->gBest.fitness)) {
                this->_failTimes=0;
                this->gBest=*curGBest;
            }
            else {
                this->_failTimes++;
            }
    }

    virtual void updatePopulation() {
        for(Particle_t & i : this->_population) {
            for(size_t idx=0;idx<Base_t::dimensions();idx++) {
                i.velocity[idx]=
                        this->_option.inertiaFactor*i.velocity[idx]
                        +this->_option.learnFactorP*randD()*(i.pBest.position[idx]-i.position[idx])
                        +this->_option.learnFactorG*randD()*(this->gBest.position[idx]-i.position[idx]);
                if(std::abs(i.velocity[idx])>this->_velocityMax[idx]) {
                    i.velocity[idx]=sign(i.velocity[idx])*this->_velocityMax[idx];
                }
                i.position[idx]+=i.velocity[idx];
                i.position[idx]=std::max(i.position[idx],this->_posMin[idx]);
                i.position[idx]=std::min(i.position[idx],this->_posMax[idx]);
            }
        }
    }

    static_assert(!(std::is_scalar<Var_t>::value),"Var_t should be a non-scalar type");

};

///Convenient typedef for stdArray (fix-sized and Runtime sized)
template<size_t DIM,
        FitnessOption FitnessOpt,
         RecordOption RecordOpt,
         class Arg_t=void,
         typename internal::PSOParameterPack<stdVecD_t<DIM>,double,Arg_t>::iFun_t _iFun_=nullptr,
         typename internal::PSOParameterPack<stdVecD_t<DIM>,double,Arg_t>::fFun_t _fFun_=nullptr>
using PSO_std = PSO<stdVecD_t<DIM>,DIM,FitnessOpt,RecordOpt,Arg_t,_iFun_,_fFun_>;


#ifdef EIGEN_CORE_H

template<size_t DIM,FitnessOption FitnessOpt,RecordOption RecordOpt,class Arg_t=void,
         typename internal::PSOParameterPack<EigenVecD_t<DIM>,double,Arg_t>::iFun_t _iFun_=nullptr,
         typename internal::PSOParameterPack<EigenVecD_t<DIM>,double,Arg_t>::fFun_t _fFun_=nullptr>
using PSO_Eigen = PSO<EigenVecD_t<DIM>,DIM,FitnessOpt,RecordOpt,Arg_t,_iFun_,_fFun_>;


///Partial specilization for PSO using Eigen's fix-sized Array
template<
        size_t DIM,
         FitnessOption FitnessOpt,
         RecordOption RecordOpt,
         class Arg_t,
        typename internal::PSOParameterPack<EigenVecD_t<DIM>,double,Arg_t>::iFun_t _iFun_,
        typename internal::PSOParameterPack<EigenVecD_t<DIM>,double,Arg_t>::fFun_t _fFun_>
class PSO<EigenVecD_t<DIM>,DIM,FitnessOpt,RecordOpt,Arg_t,_iFun_,_fFun_>
///Partial specilization for PSO using Eigen's fix-sized Array
    : public internal::PSOBase<EigenVecD_t<DIM>,DIM,double,RecordOpt,Arg_t,_iFun_,_fFun_>
{
public:
    using Base_t = internal::PSOBase<EigenVecD_t<DIM>,DIM,double,RecordOpt,Arg_t,_iFun_,_fFun_>;
    Heu_MAKE_PSOABSTRACT_TYPES
    using Var_t = EigenVecD_t<DIM>;

    static const DoubleVectorOption Flag = DoubleVectorOption::Eigen;

    virtual void setPVRange(double pMin,double pMax,double vMax) {
        this->_posMin.setConstant(this->dimensions(),1,pMin);
        this->_posMax.setConstant(this->dimensions(),1,pMax);
        this->_velocityMax.setConstant(this->dimensions(),1,vMax);
    }

    virtual double bestFitness() const {
        return this->gBest.fitness;
    }

    ///Candidate function for initialization
    static void default_iFun(Var_t *x,Var_t *v,
    const Var_t * xMin,const Var_t * xMax,
    const Var_t *,const Args_t*) {
        x->setRandom(xMin->size(),1);
        (*x)*=(*xMax-*xMin)/2;
        (*x)+=(*xMin+*xMax)/2;
        v->setZero(xMin->size(),1);
    }

protected:
    static bool isBetterThan(double a,double b) {
        if(FitnessOpt==FitnessOption::FITNESS_GREATER_BETTER) {
            return a>b;
        } else {
            return a<b;
        }
    }

    virtual void updatePGBest() {
        Point_t * curGBest=&this->_population.front().pBest;

        for(Particle_t & i : this->_population) {
            if(isBetterThan(i.fitness,i.pBest.fitness)) {
                i.pBest=i;
            }

            if(isBetterThan(i.pBest.fitness,curGBest->fitness)) {
                curGBest=&i.pBest;
            }
        }

            if(isBetterThan(curGBest->fitness,this->gBest.fitness)) {
                this->_failTimes=0;
                this->gBest=*curGBest;
            }
            else {
                this->_failTimes++;
            }
    }

    virtual void updatePopulation() {
        for(Particle_t & i : this->_population) {
            i.velocity=this->_option.inertiaFactor*i.velocity
                        +this->_option.learnFactorP*internal::randD()*(i.pBest.position-i.position)
                        +this->_option.learnFactorG*internal::randD()*(this->gBest.position-i.position);

            i.velocity=i.velocity.min(this->_velocityMax);
            i.velocity=i.velocity.max(-this->_velocityMax);

            i.position+=i.velocity;

            i.position=i.position.min(this->_posMax);
            i.position=i.position.max(this->_posMin);

        }
    }

};
#endif // EIGEN_CORE_H

}

#endif // PSO_HPP
