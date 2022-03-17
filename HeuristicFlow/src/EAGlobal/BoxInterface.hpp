// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BOXINTERFACE_HPP
#define BOXINTERFACE_HPP

#include "BoxRTCTDims.hpp"
#include "BoxReal.h"

namespace Heu
{

/**
 * Fixed(N) dim real(d) square(S) box
 */
template<size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         size_t RangeType=Runtime,
         internal::TemplateVal_t<double> MinCT=encode<0,1>::code,
         internal::TemplateVal_t<double> MaxCT=encode<0,1>::code,
         internal::TemplateVal_t<double> LRCT=encode<0,1>::code>
using BoxNdS = internal::RealBox<double,Dim,DVO,BoxShape::SQUARE_BOX,
    RangeType,MinCT,MaxCT,LRCT>;

/**
 * Runtime(X) dim real(d) square(S) box
 */
template<DoubleVectorOption DVO=DoubleVectorOption::Std,
         size_t RangeType=Runtime,
         internal::TemplateVal_t<double> MinCT=encode<0,1>::code,
         internal::TemplateVal_t<double> MaxCT=encode<0,1>::code,
         internal::TemplateVal_t<double> LRCT=encode<0,1>::code>
using BoxXdS = internal::RealBox<double,Runtime,DVO,BoxShape::SQUARE_BOX,
    RangeType,MinCT,MaxCT,LRCT>;

/**
 * Fixed(N) dim real(d) nonsquare(N) box
 */
template<size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxNdN = internal::RealBox<double,Dim,DVO,BoxShape::RECTANGLE_BOX>;

/**
 * Runtime(X) dim real(d) nonsquare(N) box
 */
template<DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxXdN = internal::RealBox<double,Runtime,DVO,BoxShape::RECTANGLE_BOX>;

/**
 * Binary box
 */
template<size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
class BooleanBox : public internal::BoxDims<bool,Dim,DVO,BoxShape::SQUARE_BOX,1,false,true>
{
public:
    static const constexpr EncodeType Encoding=EncodeType::Binary;
};

/**
 * Binary box with fixed dim
 */
template<size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxNb = typename std::enable_if<Dim!=Runtime,BooleanBox<Dim,DVO>>::type;

/**
 * Binary box with dynamic dims
 */
template<DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxXb = BooleanBox<Runtime,DVO>;


namespace internal
{
/**
 * @brief Symbolic box constraint
 */
template<typename Scalar_t,size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         BoxShape BS=BoxShape::SQUARE_BOX,size_t RangeType=Runtime,
         Scalar_t MinCT=0,Scalar_t MaxCT=1>
class SymbolBox : public internal::BoxDims<Scalar_t,Dim,DVO,BS,RangeType,MinCT,MaxCT>
{
private:
    static_assert(std::is_integral<Scalar_t>::value,"Symbol box requires integer Scalar_t");
public:
    static const constexpr EncodeType Encoding=EncodeType::Symbolic;
};

}   //  internal

/**
 * @brief Square symbolic box with fixed dim
 */
template<typename Scalar_t,size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         size_t RangeType=Runtime,
         Scalar_t MinCT=0,Scalar_t MaxCT=1>
using BoxNsS = typename std::enable_if<Dim!=Runtime,
    internal::SymbolBox<Scalar_t,Dim,DVO,BoxShape::SQUARE_BOX,
        RangeType,MinCT,MaxCT>>::type;


/**
 * @brief Square symbolic box with runtime dim
 */
template<typename Scalar_t,DoubleVectorOption DVO=DoubleVectorOption::Std,
         size_t RangeType=Runtime,
         Scalar_t MinCT=0,Scalar_t MaxCT=1>
using BoxXsS = internal::SymbolBox<Scalar_t,Runtime,DVO,
    BoxShape::SQUARE_BOX,RangeType,MinCT,MaxCT>;


/**
 * @brief Non-square symbolic box with fixed dim
 */
template<typename Scalar_t,size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxNsN = typename std::enable_if<Dim!=Runtime,
    internal::SymbolBox<Scalar_t,Dim,DVO,BoxShape::RECTANGLE_BOX>>::type;


/**
 * @brief Non-square symbolic box with runtime dim
 */
template<typename Scalar_t,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxXsN = internal::SymbolBox<Scalar_t,Runtime,DVO,BoxShape::RECTANGLE_BOX>;

/*
template<typename Scalar_t,size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         BoxShape BS=BoxShape::SQUARE_BOX,size_t RangeType=Runtime,
         Scalar_t MinCT=0,Scalar_t MaxCT=1>
class IntegerBox : public BoxDims<Scalar_t,Dim,DVO,BS,RangeType,MinCT,MaxCT>
{
private:
    static_assert(std::is_integral<Scalar_t>::value,"Use integer as Scalar_t");
public:
    static const constexpr EncodeType Encoding=EncodeType::Integer;
};

template<size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         size_t RangeType=Runtime,
         int MinCT=0,int MaxCT=1>
using BoxNiS = typename std::enable_if<Dim!=Runtime,
    IntegerBox<int,Dim,DVO,BoxShape::SQUARE_BOX,RangeType,MinCT,MaxCT>>::type;


template<DoubleVectorOption DVO=DoubleVectorOption::Std,
         size_t RangeType=Runtime,
         int MinCT=0,int MaxCT=1>
using BoxXiS = IntegerBox<int,Runtime,DVO,BoxShape::SQUARE_BOX,
    RangeType,MinCT,MaxCT>;


template<size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxNiN = typename std::enable_if<Dim!=Runtime,
    IntegerBox<int,Dim,DVO,BoxShape::RECTANGLE_BOX>>::type;



template<DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxXiN = IntegerBox<int,Runtime,DVO,BoxShape::RECTANGLE_BOX>;

*/


}   //  Heu

#endif // BOXINTERFACE_HPP
