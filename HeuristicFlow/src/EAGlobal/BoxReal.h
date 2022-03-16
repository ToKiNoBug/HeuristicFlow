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

#ifndef BOXREAL_H
#define BOXREAL_H

#include "BoxRTCTDims.hpp"

namespace Heu
{
/**
 * @brief Real boxes of different types
 */
template<typename Scalar_t,size_t Dim,
         DoubleVectorOption DVO,BoxShape BS,
         size_t RangeType,
         TemplateVal_t<Scalar_t> MinCT,
         TemplateVal_t<Scalar_t> MaxCT>
class RealBoxBase : public BoxDims<Scalar_t,Dim,DVO,BS,RangeType,MinCT,MaxCT>
{
private:
    static_assert(std::is_floating_point<Scalar_t>::value,
        "Scalar_t must be a floating point number");
public:
    static const constexpr EncodeType Encoding=EncodeType::Real;
};


template<typename Var_t>
struct learnRateBody
{
protected:
    Var_t _learnRate;
public:
    inline Var_t & learnRate() {
        return _learnRate;
    }

    inline const Var_t & learnRate() const {
        return _learnRate;
    }

    inline void setLearnRate(const Var_t & v) {
        _learnRate=v;
    }
};

/**
 * @brief Compile-time ranged box with learning rate
 */
template<typename Scalar_t,size_t Dim,
         DoubleVectorOption DVO,BoxShape BS,
         size_t RangeType=Runtime,
         TemplateVal_t<Scalar_t> MinCT=TemplateVal_t<Scalar_t>(1),
         TemplateVal_t<Scalar_t> MaxCT=TemplateVal_t<Scalar_t>(1),
         TemplateVal_t<Scalar_t> LearnRateCT=TemplateVal_t<Scalar_t>(1)>
class RealBox :
        public RealBoxBase<Scalar_t,Dim,DVO,BS,RangeType,MinCT,MaxCT>
{
private:
    static_assert(RangeType!=Runtime,"Wrong specialization of RealBox");
    static const constexpr Scalar_t learnRateCT=decode<LearnRateCT>::real;
public:
    inline constexpr Scalar_t learnRate() const {
        return learnRateCT;
    }
};

/**
 * @brief Runtime ranged box with learning rate
 */
template<typename Scalar_t,size_t Dim,
         DoubleVectorOption DVO,BoxShape BS,
         TemplateVal_t<Scalar_t> MinCT,
         TemplateVal_t<Scalar_t> MaxCT,
         TemplateVal_t<Scalar_t> LearnRateCT>
class RealBox<Scalar_t,Dim,DVO,BS,Runtime,MinCT,MaxCT,LearnRateCT>
        :
        public RealBoxBase<Scalar_t,Dim,DVO,BS,Runtime,MinCT,MaxCT> ,
        public learnRateBody<typename std::conditional<BS==BoxShape::RECTANGLE_BOX,
            Container<Scalar_t,Dim,DVO>,Scalar_t>::type>
{
};

}   //  namespace

#endif // BOXREAL_H
