// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BOXREAL_H
#define BOXREAL_H

#include "BoxRTCTDims.hpp"

namespace Eigen
{

namespace internal
{

/**
 * @brief Real boxes of different types
 */
template<typename Scalar_t,int Dim,
         DoubleVectorOption DVO,BoxShape BS,
         bool isFixedRange,
         TemplateVal_t<Scalar_t> MinCT,
         TemplateVal_t<Scalar_t> MaxCT>
class RealBoxBase : public BoxDims<Scalar_t,Dim,DVO,BS,isFixedRange,MinCT,MaxCT>
{
private:
    static_assert(std::is_floating_point<Scalar_t>::value,
        "Scalar_t must be a floating point number");
    static_assert(Dim>0||Dim==Eigen::Dynamic,"Invalid template parameter Dim");
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
template<typename Scalar_t,int Dim,
         DoubleVectorOption DVO,BoxShape BS,
         bool isFixedRange=false,
         TemplateVal_t<Scalar_t> MinCT=TemplateVal_t<Scalar_t>(1),
         TemplateVal_t<Scalar_t> MaxCT=TemplateVal_t<Scalar_t>(1),
         TemplateVal_t<Scalar_t> LearnRateCT=TemplateVal_t<Scalar_t>(1)>
class RealBox :
        public RealBoxBase<Scalar_t,Dim,DVO,BS,isFixedRange,MinCT,MaxCT>
{
private:
    static_assert(isFixedRange,"Wrong specialization of RealBox");
    static const constexpr Scalar_t learnRateCT=DivDecode<LearnRateCT>::real;
public:
    inline constexpr Scalar_t learnRate() const {
        return learnRateCT;
    }
};

/**
 * @brief runtime ranged box with learning rate
 */
template<typename Scalar_t,
         int Dim,
         DoubleVectorOption DVO,BoxShape BS,
         TemplateVal_t<Scalar_t> MinCT,
         TemplateVal_t<Scalar_t> MaxCT,
         TemplateVal_t<Scalar_t> LearnRateCT>
class RealBox<Scalar_t,Dim,DVO,BS,false,MinCT,MaxCT,LearnRateCT>
        :
        public RealBoxBase<Scalar_t,Dim,DVO,BS,false,MinCT,MaxCT> ,
        public learnRateBody<typename std::conditional<BS==BoxShape::RECTANGLE_BOX,
            Container<Scalar_t,Dim,DVO>,Scalar_t>::type>
{
};

}   //  internal

}   //  namespace Eigen

#endif // BOXREAL_H
