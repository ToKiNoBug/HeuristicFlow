#ifndef BOXINTERFACE_HPP
#define BOXINTERFACE_HPP

#include "BoxRTCTDims.hpp"

namespace Heu
{
/**
 * @brief Real boxes of different types
 */
template<typename Scalar_t,size_t Dim,
         DoubleVectorOption DVO,BoxShape BS,
         size_t RangeType=Runtime,
         TemplateVal_t<Scalar_t> MinCT=TemplateVal_t<Scalar_t>(1),
         TemplateVal_t<Scalar_t> MaxCT=TemplateVal_t<Scalar_t>(1)>
using RealBox = typename std::enable_if<std::is_floating_point<Scalar_t>::value,
    BoxDims<Scalar_t,Dim,DVO,BS,RangeType,MinCT,MaxCT>>::type;


/**
 * Fixed(N) dim real(d) square(S) box
 */
template<size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std,
         size_t RangeType=Runtime,
         TemplateVal_t<double> MinCT=encode<0,1>::code,
         TemplateVal_t<double> MaxCT=encode<0,1>::code>
using BoxNdS = RealBox<double,Dim,DVO,BoxShape::SQUARE_BOX,
    RangeType,MinCT,MaxCT>;

/**
 * Runtime(X) dim real(d) square(S) box
 */
template<DoubleVectorOption DVO=DoubleVectorOption::Std,
         size_t RangeType=Runtime,
         TemplateVal_t<double> MinCT=encode<0,1>::code,
         TemplateVal_t<double> MaxCT=encode<0,1>::code>
using BoxXdS = RealBox<double,Runtime,DVO,BoxShape::SQUARE_BOX,
    RangeType,MinCT,MaxCT>;

/**
 * Fixed(N) dim real(d) nonsquare(N) box
 */
template<size_t Dim,DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxNdN = RealBox<double,Dim,DVO,BoxShape::RECTANGLE_BOX>;

/**
 * Runtime(X) dim real(d) nonsquare(N) box
 */
template<DoubleVectorOption DVO=DoubleVectorOption::Std>
using BoxXdN = RealBox<double,Runtime,DVO,BoxShape::RECTANGLE_BOX>;

};

#endif // BOXINTERFACE_HPP
