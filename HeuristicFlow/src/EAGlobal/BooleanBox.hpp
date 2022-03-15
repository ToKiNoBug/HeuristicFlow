#ifndef Heu_BOOLEANBOX_HPP
#define Heu_BOOLEANBOX_HPP

#include "../Global/Enumerations.hpp"
#include "../Global/Constants.hpp"
#include "../Global/Types.hpp"

#include "BoxConstraint.hpp"

namespace Heu
{

template<size_t DIM,DoubleVectorOption DVO>
class BoxBooleanBase
{
public:
#ifdef EIGEN_CORE_H
    using Var_t =
        typename std::conditional
            <DVO==DoubleVectorOption::Eigen,
            stdContainer<bool,DIM>,
            EigenContainer<bool,DIM>>::type;
#else
    using Var_t = stdContainer<bool,DIM>;
    static_assert(DVO!=DoubleVectorOption::Eigen,"Please include libEigen.");
#endif  //  #ifdef EIGEN_CORE_H

    constexpr bool max() const {
        return true;
    }

    constexpr bool min() const {
        return false;
    }
    
    static const bool Heu_isBox=true;
    static const bool Heu_isSquareBox=true;
    static const EncodeType encodeType=BinaryEncoding;
};

template<size_t DIM,DoubleVectorOption DVO>
class BoxBoolean
    : public BoxBooleanBase<DIM,DVO>
{
public:
    using Base_t = BoxBooleanBase<DIM,DVO>;
    Heu_MAKE_BOXABSTRACT_TYPES

    inline constexpr size_t varDim() const {
        return DIM;
    }

    inline void setVarDim(size_t d) const {
        assert(d==DIM);
    }
};


template<DoubleVectorOption DVO>
class BoxBoolean<Runtime,DVO>
    : public BoxBooleanBase<Runtime,DVO>
{
public:
    using Base_t = BoxBooleanBase<Runtime,DVO>;
    Heu_MAKE_BOXABSTRACT_TYPES

    BoxBoolean() {
        _varDim=3;
    }

    inline constexpr size_t varDim() const {
        return _varDim;
    }

    inline void setVarDim(size_t d) const {
        _varDim=d;
    }
protected:
    size_t _varDim;
};

namespace HeuPrivate
{

template<typename Scalar_t,size_t DIM,BoxShape BS,DoubleVectorOption DVO>
struct box_traits
{
private:
    static const bool isReal=std::is_floating_point<Scalar_t>::value;
    static const bool isBoolean=std::is_same<Scalar_t,bool>::value;
    static const bool isSymbolic=(!isReal)&&(!isBoolean);

    static_assert(isReal||(BS==SQUARE_BOX),"Non-real box must be a square box");
public:
    using type = typename std::conditional<
        isReal,
        BoxFloat<Scalar_t,DIM,BS,DVO>,
        typename std::conditional<
            isBoolean,BoxBoolean<DIM,DVO>,
            BoxSymbolic<Scalar_t,DIM,DVO>
            >::type
        >::type;
};

template<typename Scalar_t,size_t DIM,DoubleVectorOption DVO,BoxShape BS=SQUARE_BOX>
using pri_BoxConstraint=typename box_traits<Scalar_t,DIM,BS,DVO>::type;

template<typename Scalar_t,size_t DIM,BoxShape BS,DoubleVectorOption DVO,class Args_t>
struct box_parameterPack
    : public pri_BoxConstraint<Scalar_t,DIM,DVO,BS>
{
    using Base_t = pri_BoxConstraint<Scalar_t,DIM,DVO,BS>;
    Heu_MAKE_BOXABSTRACT_TYPES
    
    Args_t & arg() {
        return _arg;
    }

    const Args_t & arg() const {
        return _arg;
    }

protected:
    Args_t _arg;
private:
    static_assert(!std::is_same<Args_t,void>::value,
        "Args_t can't be void");
};}   //  namespace HeuPrivate


template<typename Scalar_t,size_t DIM,DoubleVectorOption DVO=
#ifdef EIGEN_CORE_H
    DoubleVectorOption::Eigen
#else
    DoubleVectorOption::Std
#endif  //  EIGEN_CORE_H
    ,class Args_t=void,BoxShape BS=SQUARE_BOX>
using BoxConstraint = typename std::conditional<
    std::is_same<Args_t,void>::value,
    HeuPrivate::pri_BoxConstraint<Scalar_t,DIM,DVO,BS>,
    HeuPrivate::box_parameterPack<Scalar_t,DIM,BS,DVO,Args_t>
    >::type;

}   //  namespace Heu

#endif  //  Heu_BOOLEANBOX_HPP
