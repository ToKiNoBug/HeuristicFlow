/*
 Copyright © 2022  TokiNoBug
This file is part of Heuristic.

    Heuristic is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Heuristic is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Heuristic.  If not, see <https://www.gnu.org/licenses/>.

*/

#ifndef BOXSHAPES_HPP
#define BOXSHAPES_HPP

#include "BoxRTCTRange.hpp"

namespace Heu
{

template<typename Scalar_t,size_t Size,
         DoubleVectorOption DVO>
using NonsquareBox =
    BoxDynamicRange<Scalar_t,Size,BoxShape::RECTANGLE_BOX,DVO>;


/**
 * @brief Square box with runtime range
 */
template<typename Scalar_t,size_t Size,
         DoubleVectorOption DVO,
         size_t RangeType,TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
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
template<typename Scalar_t,size_t Size,
         DoubleVectorOption DVO,
         TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class SquareBox<Scalar_t,Size,DVO,Runtime,MinCT,MaxCT>
        : public BoxDynamicRange<Scalar_t,Size,BoxShape::SQUARE_BOX,DVO>
{
private:
    using Base_t = BoxDynamicRange<Scalar_t,Size,BoxShape::SQUARE_BOX,DVO>;
public:
    using Var_t = typename Base_t::Var_t;
};


}   //  namespace

#endif // BOXSHAPES_HPP
