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
template<int Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         bool isFixedRange=false,
         internal::TemplateVal_t<double> MinCT=DivEncode<0,1>::code,
         internal::TemplateVal_t<double> MaxCT=DivEncode<0,1>::code,
         internal::TemplateVal_t<double> LRCT=DivEncode<0,1>::code>
using BoxNdS = internal::RealBox<double,Dim,DVO,BoxShape::SQUARE_BOX,
    isFixedRange,MinCT,MaxCT,LRCT>;

/**
 * runtime(X) dim real(d) square(S) box
 */
template<DoubleVectorOption DVO=DoubleVectorOption::Std,
         bool isFixedRange=false,
         internal::TemplateVal_t<double> MinCT=DivEncode<0,1>::code,
         internal::TemplateVal_t<double> MaxCT=DivEncode<0,1>::code,
         internal::TemplateVal_t<double> LRCT=DivEncode<0,1>::code>
using BoxXdS = internal::RealBox<double,Eigen::Dynamic,DVO,BoxShape::SQUARE_BOX,
    isFixedRange,MinCT,MaxCT,LRCT>;

/**
 * Fixed(N) dim real(d) nonsquare(N) box
 */
template<int Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxNdN = internal::RealBox<double,Dim,DVO,BoxShape::RECTANGLE_BOX>;

/**
 * runtime(X) dim real(d) nonsquare(N) box
 */
template<DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxXdN = internal::RealBox<double,Eigen::Dynamic,DVO,BoxShape::RECTANGLE_BOX>;

/**
 * Binary box
 */
template<int Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
class BooleanBox : public internal::BoxDims<bool,Dim,DVO,BoxShape::SQUARE_BOX,true,0,1>
{
private:
    static_assert(Dim>0||Dim==Eigen::Dynamic,"Invalid template parameter Dim");
public:
    static const constexpr EncodeType Encoding=EncodeType::Binary;
};

/**
 * Binary box with fixed dim
 */
template<int Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxNb = typename std::enable_if<Dim!=Eigen::Dynamic,BooleanBox<Dim,DVO>>::type;

/**
 * Binary box with dynamic dims
 */
template<DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxXb = BooleanBox<Eigen::Dynamic,DVO>;


namespace internal
{
/**
 * @brief Symbolic box constraint
 */
template<typename Scalar_t,int Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         BoxShape BS=BoxShape::SQUARE_BOX,bool isFixedRange=false,
         Scalar_t MinCT=0,Scalar_t MaxCT=1>
class SymbolBox : public internal::BoxDims<Scalar_t,Dim,DVO,BS,isFixedRange,MinCT,MaxCT>
{
private:
    static_assert(std::is_integral<Scalar_t>::value,"Symbol box requires integer Scalar_t");
    static_assert(Dim>0||Dim==Eigen::Dynamic,"Invalid template parameter Dim");
public:
    static const constexpr EncodeType Encoding=EncodeType::Symbolic;
};

}   //  internal

/**
 * @brief Square symbolic box with fixed dim
 */
template<typename Scalar_t,int Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         bool isFixedRange=false,
         Scalar_t MinCT=0,Scalar_t MaxCT=1>
using BoxNsS = typename std::enable_if<Dim!=Eigen::Dynamic,
    internal::SymbolBox<Scalar_t,Dim,DVO,BoxShape::SQUARE_BOX,
        isFixedRange,MinCT,MaxCT>>::type;


/**
 * @brief Square symbolic box with runtime dim
 */
template<typename Scalar_t,DoubleVectorOption DVO=DoubleVectorOption::Std,
         bool isFixedRange=false,
         Scalar_t MinCT=0,Scalar_t MaxCT=1>
using BoxXsS = internal::SymbolBox<Scalar_t,Eigen::Dynamic,DVO,
    BoxShape::SQUARE_BOX,isFixedRange,MinCT,MaxCT>;


/**
 * @brief Non-square symbolic box with fixed dim
 */
template<typename Scalar_t,int Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxNsN = typename std::enable_if<Dim!=Eigen::Dynamic,
    internal::SymbolBox<Scalar_t,Dim,DVO,BoxShape::RECTANGLE_BOX>>::type;


/**
 * @brief Non-square symbolic box with runtime dim
 */
template<typename Scalar_t,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxXsN = internal::SymbolBox<Scalar_t,Eigen::Dynamic,DVO,BoxShape::RECTANGLE_BOX>;


}   //  Heu

#endif // BOXINTERFACE_HPP
