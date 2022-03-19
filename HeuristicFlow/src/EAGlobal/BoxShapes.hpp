// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BOXSHAPES_HPP
#define BOXSHAPES_HPP

#include <Eigen/Core>
#include "BoxRTCTRange.hpp"

namespace Heu
{

namespace internal
{

template<typename Scalar_t,int Size,
         DoubleVectorOption DVO>
using NonsquareBox =
    BoxDynamicRange<Scalar_t,Size,BoxShape::RECTANGLE_BOX,DVO>;


/**
 * @brief Square box with runtime range
 */
template<typename Scalar_t,int Size,
         DoubleVectorOption DVO,
         bool isFixedRange,TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class SquareBox : public BoxFixedRange<Scalar_t,MinCT,MaxCT>
{
private:
    using Base_t = BoxFixedRange<Scalar_t,MinCT,MaxCT>;
public:
    using Var_t = Container<Scalar_t,Size,DVO>;
    //using Var_t = typename Base_t::Var_t;

protected:

};

/**
 * @brief Square box with runtime ranges
 */
template<typename Scalar_t,int Size,
         DoubleVectorOption DVO,
         TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class SquareBox<Scalar_t,Size,DVO,false,MinCT,MaxCT>
        : public BoxDynamicRange<Scalar_t,Size,BoxShape::SQUARE_BOX,DVO>
{
private:
    using Base_t = BoxDynamicRange<Scalar_t,Size,BoxShape::SQUARE_BOX,DVO>;
public:
    using Var_t = typename Base_t::Var_t;
};

}   //  internal

}   //  namespace

#endif // BOXSHAPES_HPP
