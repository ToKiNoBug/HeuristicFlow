/*
 Copyright © 2022  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef PSO_HPP
#define PSO_HPP
#include "PSOOption.hpp"
#include "PSOBase.hpp"
#include <array>
#include <vector>
#include <tuple>
#include <type_traits>

namespace OptimT {

///Generalized PSO solver
template<class Var_t,
         size_t DIM,
         FitnessOption FitnessOpt,
         RecordOption RecordOpt,
         class ...Args>
class PSO : public PSOBase<Var_t,DIM,double,RecordOpt,Args...>
{
public:
    using Base_t = PSOBase<Var_t,DIM,double,RecordOpt,Args...>;
    OPTIMT_MAKE_PSOABSTRACT_TYPES

    static const DoubleVectorOption Flag =
        (std::is_same<Var_t,stdVecD_t<DIM>>::value)?
        DoubleVectorOption::Std : DoubleVectorOption::Custom;
    
    void setPVRange(double pMin,double pMax,double vMax) {
        std::cerr<<__FILE__<<" , "<<__LINE__<<std::endl;
        for(size_t i=0;i<this->dimensions();i++) {
            Base_t::_posMin[i]=pMin;
            Base_t::_posMax[i]=pMax;
            Base_t::_velocityMax[i]=vMax;
        }
    }
    
    ///function used to provide a result for recording
    virtual double bestFitness() const {
        return Base_t::gBest.fitness;
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
        Point_t * curGBest=&Base_t::_population.front().pBest;

        for(Particle_t & i : Base_t::_population) {
            if(isBetterThan(i.fitness,i.pBest.fitness)) {
                i.pBest=i;
            }

            if(isBetterThan(i.pBest.fitness,curGBest->fitness)) {
                curGBest=&i.pBest;
            }
        }

            if(isBetterThan(curGBest->fitness,Base_t::gBest.fitness)) {
                Base_t::_failTimes=0;
                Base_t::gBest=*curGBest;
            }
            else {
                Base_t::_failTimes++;
            }
    }

    virtual void updatePopulation() {
        for(Particle_t & i : Base_t::_population) {
            for(size_t idx=0;idx<Base_t::dimensions();idx++) {
                i.velocity[idx]=
                        Base_t::_option.inertiaFactor*i.velocity[idx]
                        +Base_t::_option.learnFactorP*randD()*(i.pBest.position[idx]-i.position[idx])
                        +Base_t::_option.learnFactorG*randD()*(Base_t::gBest.position[idx]-i.position[idx]);
                if(std::abs(i.velocity[idx])>Base_t::_velocityMax[idx]) {
                    i.velocity[idx]=sign(i.velocity[idx])*Base_t::_velocityMax[idx];
                }
                i.position[idx]+=i.velocity[idx];
                i.position[idx]=std::max(i.position[idx],Base_t::_posMin[idx]);
                i.position[idx]=std::min(i.position[idx],Base_t::_posMax[idx]);
            }
        }
    }

    static_assert(!(std::is_scalar<Var_t>::value),"Var_t should be a non-scalar type");

};

///Convenient typedef for stdArray (fix-sized and dynamic sized)
template<size_t DIM,
        FitnessOption FitnessOpt,
         RecordOption RecordOpt,
         class ...Args>
using PSO_std = PSO<stdVecD_t<DIM>,DIM,FitnessOpt,RecordOpt,Args...>;


#ifdef OptimT_PSO_USE_EIGEN

///Convenient typedef for Eigen's Array (fix-sized and dynamic-sized)
template<size_t DIM,FitnessOption FitnessOpt,RecordOption RecordOpt,class ...Args>
using PSO_Eigen = PSO<EigenVecD_t<DIM>,DIM,FitnessOpt,RecordOpt,Args...>;
///Convenient typedef for Eigen's Array (fix-sized and dynamic-sized)


///Partial specilization for PSO using Eigen's fix-sized Array
template<
        size_t DIM,
         FitnessOption FitnessOpt,
         RecordOption RecordOpt,
         class ...Args>
class PSO<EigenVecD_t<DIM>,DIM,FitnessOpt,RecordOpt,Args...>
///Partial specilization for PSO using Eigen's fix-sized Array
    : public PSOBase<EigenVecD_t<DIM>,DIM,double,RecordOpt,Args...>
{
public:
    using Base_t = PSOBase<EigenVecD_t<DIM>,DIM,double,RecordOpt,Args...>;
    OPTIMT_MAKE_PSOABSTRACT_TYPES
    using Var_t = EigenVecD_t<DIM>;

    static const DoubleVectorOption Flag = DoubleVectorOption::Eigen;

    virtual void setPVRange(double pMin,double pMax,double vMax) {
        this->_posMin.setConstant(this->dimensions(),1,pMin);
        this->_posMax.setConstant(this->dimensions(),1,pMax);
        this->_velocityMax.setConstant(this->dimensions(),1,vMax);
    }

    virtual double bestFitness() const {
        return Base_t::gBest.fitness;
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
        Point_t * curGBest=&Base_t::_population.front().pBest;

        for(Particle_t & i : Base_t::_population) {
            if(isBetterThan(i.fitness,i.pBest.fitness)) {
                i.pBest=i;
            }

            if(isBetterThan(i.pBest.fitness,curGBest->fitness)) {
                curGBest=&i.pBest;
            }
        }

            if(isBetterThan(curGBest->fitness,Base_t::gBest.fitness)) {
                Base_t::_failTimes=0;
                Base_t::gBest=*curGBest;
            }
            else {
                Base_t::_failTimes++;
            }
    }

    virtual void updatePopulation() {
        for(Particle_t & i : Base_t::_population) {
            i.velocity=this->_option.inertiaFactor*i.velocity
                        +this->_option.learnFactorP*randD()*(i.pBest.position-i.position)
                        +this->_option.learnFactorG*randD()*(this->gBest.position-i.position);

            i.velocity=i.velocity.min(Base_t::_velocityMax);
            i.velocity=i.velocity.max(-Base_t::_velocityMax);

            i.position+=i.velocity;

            i.position=i.position.min(Base_t::_posMax);
            i.position=i.position.max(Base_t::_posMin);

        }
    }

};
#endif // OptimT_PSO_USE_EIGEN

}

#endif // PSO_HPP
