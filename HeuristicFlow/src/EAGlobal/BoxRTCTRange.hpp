// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BOX_RTCTRange_HPP
#define BOX_RTCTRange_HPP
#include <stdint.h>
#include <assert.h>
#include "../../Global"

/**
 * This file has implementation for different min() and max()
 * and incomplete setMin() and setMax() implementation
 */

namespace Eigen
{

namespace internal
{
/**
 * @brief Non-square box with runtime range.
 */
template<typename Scalar_t,int Size,BoxShape BS,DoubleVectorOption DVO>
class BoxDynamicRange
{
public:
    using Var_t = Container<Scalar_t,Size,DVO>;
    static const constexpr bool isBox=true;
    static const constexpr char Flag[]="Non-square box with runtime range.";
    static const constexpr BoxShape Shape = BoxShape::RECTANGLE_BOX;

    inline Var_t & min() {
        return _minV;
    }

    inline Var_t & max() {
        return _maxV;
    }

    inline const Var_t & min() const {
        return _minV;
    }

    inline const Var_t & max() const {
        return _maxV;
    }

    inline void setMin(const Var_t _m) {
        _minV=_m;
    }

    inline void setMax(const Var_t _m) {
        _maxV=_m;
    }
protected:

    Var_t _minV;
    Var_t _maxV;
private:
    static_assert (BS!=BoxShape::SQUARE_BOX,
        "Non square box specialization used for square box");
};


/**
 * @brief Square box with runtime range
 */
template<typename Scalar_t,int Size,DoubleVectorOption DVO>
class BoxDynamicRange<Scalar_t,Size,BoxShape::SQUARE_BOX,DVO>
{
public:
    static const constexpr bool isBox=true;
    static const constexpr char Flag[]="Square box with runtime range.";
    static const constexpr BoxShape Shape = BoxShape::SQUARE_BOX;
    using Var_t = Container<Scalar_t,Size,DVO>;

    inline Scalar_t min() const {
        return _minS;
    }

    inline Scalar_t max() const {
        return _maxS;
    }

    inline Scalar_t & min() {
        return _minS;
    }

    inline Scalar_t & max() {
        return _maxS;
    }

    inline void setMin(Scalar_t s) {
        _minS=s;
    }

    inline void setMax(Scalar_t s) {
        _maxS=s;
    }

protected:

    Scalar_t _minS;
    Scalar_t _maxS;
};


template<typename Scalar_t>
using TemplateVal_t = typename std::conditional<
    std::is_floating_point<Scalar_t>::value,DivCode,Scalar_t>::type;

/**
 * @brief Square box with compile-time range
 */
template<typename Scalar_t,
         TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class BoxFixedRange
{
public:
    static const constexpr bool isBox=true;
    static const constexpr char Flag[]="Square box with compile-time range";
    static const constexpr BoxShape Shape = BoxShape::SQUARE_BOX;
    inline constexpr Scalar_t min() const {
        return Decoder<isFloatPoint>::minCT;
    }

    inline constexpr Scalar_t max() const {
        return Decoder<isFloatPoint>::maxCT;
    }
protected:

private:
    static const constexpr bool isFloatPoint=std::is_floating_point<Scalar_t>::value;

    template<bool isFl,typename unused=void>
    struct Decoder
    {
        static const constexpr Scalar_t minCT=MinCT;
        static const constexpr Scalar_t maxCT=MaxCT;
        static_assert(isFl==isFloatPoint,"Wrong specialization");
    };

    template<typename unused>
    struct Decoder<true,unused>
    {
        static const constexpr Scalar_t minCT=DivDecode<MinCT>::real;
        static const constexpr Scalar_t maxCT=DivDecode<MaxCT>::real;
        static_assert(isFloatPoint,"Wrong specialization");
    };
};

}   //  internal

}   //  namespace Eigen

#endif // BOX_RTCTRange_HPP
