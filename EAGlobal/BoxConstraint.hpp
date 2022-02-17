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

#ifndef Heu_BOXCONSTRAINT_HPP
#define Heu_BOXCONSTRAINT_HPP

#include "../Global/Enumerations.hpp"
#include "../Global/Constants.hpp"
#include "../Global/Types.hpp"
#include "../SimpleMatrix/MatrixDynamicSize.hpp"
#include "../SimpleMatrix/MatrixFixedSize.hpp"

#include "../Global/Globals.hpp"
#include <algorithm>

namespace Heu
{

/**
 * @brief Fundamental types defination and static assertions
 * 
 * @tparam Scalar_t Type of element
 * @tparam DIM Size of Var_t
 * @tparam DVO Type of container
 */
template<typename Scalar_t,size_t DIM,DoubleVectorOption DVO>
class BoxAbstract
{
public:
#ifdef EIGEN_CORE_H
    using Var_t =
        typename std::conditional
            <DVO==DoubleVectorOption::Eigen,
            stdContainer<Scalar_t,DIM>,
            EigenContainer<Scalar_t,DIM>>::type;
#else
    using Var_t = stdContainer<Scalar_t,DIM>;
    static_assert(DVO!=DoubleVectorOption::Eigen,"Please include libEigen.");
#endif  //  #ifdef EIGEN_CORE_H

private:
    static_assert(DVO!=DoubleVectorOption::Custom,
        "Box constraint for custom array types isn't provided");
    static_assert(DIM!=1,
        "Why use evolutionary algorithms for 1-DIM functions?");
};


#define Heu_MAKE_BOXABSTRACT_TYPES \
using Var_t = typename Base_t::Var_t;


/**
 * @brief For non-square box
 * 
 * @tparam BS Box type
 */
template<typename Scalar_t,size_t DIM,BoxShape BS,DoubleVectorOption DVO>
class BoxBase
    : public BoxAbstract<Scalar_t,DIM,DVO>
{
public:
    BoxBase() {}
    ~BoxBase() {}

    using Base_t = BoxAbstract<Scalar_t,DIM,DVO>;
    Heu_MAKE_BOXABSTRACT_TYPES


    inline void setMin(const Var_t & min) {
        _min=min;
    }

    inline void setMax(const Var_t & max) {
        _max=max;
    }

    inline void setMin(const Scalar_t s) {
        for(auto & i : _min) {
            i=s;
        }
    }

    inline void setMax(const Scalar_t s) {
        for(auto & i : _max) {
            i=s;
        }
    }

    inline const Var_t & min() const {
        return _min;
    }

    inline const Var_t & max() const {
        return _max;
    }

    static const bool Heu_isInheritedFromBoxBase=true;
    static const bool Heu_isSquareBox=false;

protected:
    Var_t _min;
    Var_t _max;
};



/**
 * @brief For square box
 * 
 */
template<typename Scalar_t,size_t DIM,DoubleVectorOption DVO>
class BoxBase<Scalar_t,DIM,BoxShape::SQUARE_BOX,DVO>
    : public BoxAbstract<Scalar_t,DIM,DVO>
{
public:
    BoxBase() {
        _min_s=0;
        _max_s=1;
    }

    ~BoxBase() {}

    using Base_t = BoxAbstract<Scalar_t,DIM,DVO>;
    Heu_MAKE_BOXABSTRACT_TYPES

    inline void setMin(const Scalar_t s) {
        _min_s=s;
    }

    inline void setMax(const Scalar_t s) {
        _max_s=s;
    }

    inline Scalar_t min() const {
        return _min_s;
    }

    inline Scalar_t max() const {
        return _max_s;
    }

    static const bool Heu_isInheritedFromBoxBase=true;
    static const bool Heu_isSquareBox=true;
protected:
    Scalar_t _min_s;
    Scalar_t _max_s;
};



/**
 * @brief For fixed dimensional box
 */
template<typename Scalar_t,size_t DIM,BoxShape BS=SQUARE_BOX,DoubleVectorOption DVO=Std>
class Box
    : public BoxBase<Scalar_t,DIM,BS,DVO>
{
public:
    using Base_t = BoxBase<Scalar_t,DIM,BS,DVO>;
    Heu_MAKE_BOXABSTRACT_TYPES

    inline constexpr size_t varDim() const {
        return DIM;
    }

protected:
    static_assert(DIM!=0,"class Box's partial specialization for dynamic-sized Var_t isn't activated");
};



/**
 * @brief For runtime dimensional box
 */
template<typename Scalar_t,BoxShape BS,DoubleVectorOption DVO>
class Box<Scalar_t,Dynamic,BS,DVO>
    : public BoxBase<Scalar_t,Dynamic,BS,DVO>
{
public:
    Box() {
        setVarDim(3);
    }

    ~Box() {

    }

    using Base_t = BoxBase<Scalar_t,Dynamic,BS,DVO>;
    Heu_MAKE_BOXABSTRACT_TYPES

    inline size_t varDim() const {
        return _varDim;
    }

    inline void setVarDim(size_t v) {
        _varDim=v;
        if constexpr (BS==RECTANGLE_BOX) {
            if constexpr (DVO==Std) {
                this->_max.resize(_varDim);
                this->_min.resize(_varDim);
            }
            else {
                this->_max.resize(_varDim,1);
                this->_min.resize(_varDim,1);
            }
        }        
    }

protected:
    size_t _varDim;
};



/**
 * @brief For continous space
 */
template<typename Scalar_t,size_t DIM,BoxShape BS=SQUARE_BOX,DoubleVectorOption DVO=Std>
class BoxFloat
    : public Box<Scalar_t,DIM,BS,DVO>
{
public:
    BoxFloat() {
        _learningRate=0.1;
    }

    ~BoxFloat() {

    }

    using Base_t = Box<Scalar_t,DIM,BS,DVO>;
    Heu_MAKE_BOXABSTRACT_TYPES

    inline Scalar_t learningRate() const {
        return _learningRate;
    }

    inline void setLearningRate(Scalar_t lr) {
        assert(lr>=0);
        assert(lr<=1);
        _learningRate=lr;
    }

    static const bool Heu_isBoxFloat=true;
protected:
    Scalar_t _learningRate;

private:
    static_assert(std::is_floating_point<Scalar_t>::value,"Scalar_t isn't a floating point type");
};

template<typename Scalar_t,size_t DIM,DoubleVectorOption DVO=Std>
class BoxSymbolic
    : public Box<Scalar_t,DIM,SQUARE_BOX,DVO>
{
public:
    using Base_t = Box<Scalar_t,DIM,SQUARE_BOX,DVO>;
    Heu_MAKE_BOXABSTRACT_TYPES

    using EditMat_t = typename std::conditional<
        DIM==Dynamic,
        MatrixDynamicSize<Scalar_t>,
        MatrixFixedSize<Scalar_t,DIM,DIM-1>>::type;

    template<bool randShuffle=false>
    void makeEidtMat() {
        const size_t span=this->_max_s-this->_min_s+1;
        _editMat.resize(span,span-1);
        std::vector<Scalar_t> temp;
        temp.resize(span-1);
        for(Scalar_t r=this->_min_s;r<=this->_max_s;r++) {
            size_t idx=0;
            for(Scalar_t c=this->_min_s;c<=this->_max_s;c++) {
                if(r==c) continue;
                temp[idx++]=c;
            }
            
            if(randShuffle) {
                std::shuffle(temp.begin(),temp.end(),global_mt19937);
            }

            for(size_t c=0;c<span-1;c++) {
                _editMat(r-this->_min_s,c)=temp[c];
            }
        }
    }

    const EditMat_t & editMat() const {
        return _editMat;
    }

protected:
    EditMat_t _editMat;
};

}   //  namespace Heu

#endif  //  Heu_BOXCONSTRAINT_HPP
