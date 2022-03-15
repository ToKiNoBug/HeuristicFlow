#ifndef BOXRTCTDIMS_HPP
#define BOXRTCTDIMS_HPP

#include "BoxShapes.hpp"

namespace Heu
{

/**
 * @brief Square box with fixed dims
 */
template<typename Scalar_t,size_t Dim,
         DoubleVectorOption DVO,BoxShape BS,
         size_t RangeType,TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class BoxCTDim : public SquareBox<Scalar_t,Dim,DVO,RangeType,MinCT,MaxCT>
{
private:
    static_assert(Dim!=Runtime,"BoxFixedSize used for runtime sized box");
    static_assert(BS==BoxShape::SQUARE_BOX,"NonSquarebox used for squarebox");
    using Base_t = SquareBox<Scalar_t,Dim,DVO,RangeType,MinCT,MaxCT>;
public:
    using Var_t = typename Base_t::Var_t;

    inline void setMin(Scalar_t s) {
        Base_t::setMin(s);
    }

    inline void setMax(Scalar_t s) {
        Base_t::setMax(s);
    }

    inline constexpr size_t dimensions() const {
        return Dim;
    }
};

/**
 * @brief NonSquarebox with fixed dims
 */
template<typename Scalar_t,size_t Dim,
         DoubleVectorOption DVO,
         size_t RangeType,TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class BoxCTDim<Scalar_t,Dim,DVO,BoxShape::RECTANGLE_BOX,
        RangeType,MinCT,MaxCT> :
        public NonsquareBox<Scalar_t,Dim,DVO>
{
private:
    static_assert(RangeType==Runtime,
        "Compile-time box range for non-square box is not supported");
    static_assert(Dim!=Runtime,"BoxFixedSize used for runtime sized box");
    using Base_t = NonsquareBox<Scalar_t,Dim,DVO>;
public:
    using Var_t = typename Base_t::Var_t;

    inline void setMin(const Var_t & v) {
        Base_t::setMin(v);
    }

    inline void setMax(const Var_t & v) {
        Base_t::setMax(v);
    }

    inline constexpr size_t dimensions() const {
        return Dim;
    }
};


/**
 * @brief Square box with runtime dims
 */
template<typename Scalar_t,
         DoubleVectorOption DVO,BoxShape BS,
         size_t RangeType,TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class BoxRTDim : public SquareBox<Scalar_t,Runtime,DVO,RangeType,MinCT,MaxCT>
{
private:
    static_assert(BS==BoxShape::SQUARE_BOX,"NonSquarebox used for squarebox");
    using Base_t = SquareBox<Scalar_t,Runtime,DVO,RangeType,MinCT,MaxCT>;
protected:
    size_t _dims;
public:
    BoxRTDim() {
        _dims=0;
    }

    inline void setMin(Scalar_t s) {
        Base_t::setMin(s);
    }

    inline void setMax(Scalar_t s) {
        Base_t::setMax(s);
    }

    inline void setDimensions(size_t d) {
        assert(d!=Runtime);
        _dims=d;
    }

    inline size_t dimensions() const {
        return _dims;
    }

};

/**
 * @brief Non-square box with dynamic dims
 */
template<typename Scalar_t,
         DoubleVectorOption DVO,
         size_t RangeType,TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
class BoxRTDim<Scalar_t,DVO,BoxShape::RECTANGLE_BOX,RangeType,MinCT,MaxCT>
        : public NonsquareBox<Scalar_t,Runtime,DVO>
{
private:
    static_assert(RangeType==Runtime,
        "Compile-time box range for non-square box is not supported");
    using Base_t = NonsquareBox<Scalar_t,Runtime,DVO>;
public:
    using Var_t = typename Base_t ::Var_t;

    inline void setMin(const Var_t & v) {
        Base_t::setMin(v);
    }

    inline void setMax(const Var_t & v) {
        Base_t::setMax(v);
    }

    inline size_t dimensions() const {
        assert(this->_minV.size()==this->_maxV.size());
        return this->_minV.size();
    }
};



template<typename Scalar_t,size_t Dim,
         DoubleVectorOption DVO,BoxShape BS,
         size_t RangeType,TemplateVal_t<Scalar_t> MinCT,TemplateVal_t<Scalar_t> MaxCT>
using BoxDims = typename std::conditional<Dim==Runtime,
    BoxRTDim<Scalar_t,DVO,BS,RangeType,MinCT,MaxCT>,
    BoxCTDim<Scalar_t,Dim,DVO,BS,RangeType,MinCT,MaxCT>>::type;

}   //  namespace

#endif // BOXRTCTDIMS_HPP
