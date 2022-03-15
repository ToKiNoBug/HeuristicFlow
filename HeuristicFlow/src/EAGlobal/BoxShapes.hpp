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
class SquareBox : public BoxDynamicRange<Scalar_t,Size,BoxShape::SQUARE_BOX,DVO>
{
private:
    using Base_t = BoxDynamicRange<Scalar_t,Size,BoxShape::SQUARE_BOX,DVO>;
public:

protected:

};

/**
 * @brief Square box with compile-time ranges
 */
template<typename Scalar_t,size_t Size,
         DoubleVectorOption DVO,
         TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class SquareBox<Scalar_t,Size,DVO,Runtime,MinCT,MaxCT>
        : public BoxFixedRange<Scalar_t,MinCT,MaxCT>
{
public:
    using Var_t = Container<Scalar_t,Size,DVO>;
private:
    using Base_t = BoxFixedRange<Scalar_t,MinCT,MaxCT>;
};


/*
template<typename Scalar_t,size_t Size,
         DoubleVectorOption DVO,BoxShape BS,
         size_t RangeType,TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
using ShapedBox = typename std::conditional<BS==BoxShape::SQUARE_BOX,
    SquareBox<Scalar_t,Size,DVO,RangeType,MinCT,MaxCT>,
    NonsquareBox<Scalar_t,Size,DVO,RangeType,MinCT,MaxCT>>::type;

*/
}   //  namespace

#endif // BOXSHAPES_HPP
