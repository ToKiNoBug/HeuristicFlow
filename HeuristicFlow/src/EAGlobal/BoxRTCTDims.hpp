// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_HEU_BOXRTCTDIMS_HPP
#define EIGEN_HEU_BOXRTCTDIMS_HPP

#include "BoxShapes.hpp"

namespace Eigen
{

namespace internal
{

/**
 * @brief Square box with fixed dims
 */
template<typename Scalar_t,
  int Dim,
  DoubleVectorOption DVO,
  BoxShape BS,
  bool isFixedRange,
  TemplateVal_t<Scalar_t> MinCT,
  TemplateVal_t<Scalar_t> MaxCT>
class BoxCTDim : public SquareBox<Scalar_t,Dim,DVO,isFixedRange,MinCT,MaxCT>
{
  private:
    static_assert(Dim!=Eigen::Dynamic,"BoxFixedSize used for runtime sized box");
    static_assert(BS==BoxShape::SQUARE_BOX,"NonSquarebox used for squarebox");

    using Base_t = SquareBox<Scalar_t,Dim,DVO,isFixedRange,MinCT,MaxCT>;

  public:
    using Var_t = typename Base_t::Var_t;

    inline constexpr int dimensions() const
    {
      return Dim;
    }
};

/**
 * @brief NonSquarebox with fixed dims
 */
template<typename Scalar_t,
  int Dim,DoubleVectorOption DVO,
  bool isFixedRange,
  TemplateVal_t<Scalar_t> MinCT,
  TemplateVal_t<Scalar_t> MaxCT>
class BoxCTDim<Scalar_t,Dim,DVO,BoxShape::RECTANGLE_BOX,isFixedRange,MinCT,MaxCT>
  : public NonsquareBox<Scalar_t,Dim,DVO>
{
  private:
    static_assert(isFixedRange==false,
    "Compile-time box range for non-square box is not supported");
    static_assert(Dim!=Eigen::Dynamic,"BoxFixedSize used for runtime sized box");
    using Base_t = NonsquareBox<Scalar_t,Dim,DVO>;
  public:
    using Var_t = typename Base_t::Var_t;

    inline constexpr int dimensions() const
    {
      return Dim;
    }

};


/**
 * @brief Square box with runtime dims
 */
template<typename Scalar_t,
  DoubleVectorOption DVO,
  BoxShape BS,
  bool isFixedRange,
  TemplateVal_t<Scalar_t> MinCT,
  TemplateVal_t<Scalar_t> MaxCT>
class BoxRTDim : public SquareBox<Scalar_t,Eigen::Dynamic,DVO,isFixedRange,MinCT,MaxCT>
{
  private:
    static_assert(BS==BoxShape::SQUARE_BOX,"NonSquarebox used for squarebox");
    using Base_t = SquareBox<Scalar_t,Eigen::Dynamic,DVO,isFixedRange,MinCT,MaxCT>;
  protected:
    int _dims;
  public:
    BoxRTDim()
    {
      _dims=0;
    }

    inline void setMin(Scalar_t s)
    {
      Base_t::setMin(s);
    }

    inline void setMax(Scalar_t s)
    {
      Base_t::setMax(s);
    }

    inline void setDimensions(int d)
    {
      assert(d!=Eigen::Dynamic);
      _dims=d;
    }

    inline int dimensions() const
    {
      return _dims;
    }

};

/**
 * @brief Non-square box with dynamic dims
 */
template<typename Scalar_t,
  DoubleVectorOption DVO,
  bool isFixedRange,
  TemplateVal_t<Scalar_t> MinCT,
  TemplateVal_t<Scalar_t> MaxCT>
class BoxRTDim<Scalar_t,DVO,BoxShape::RECTANGLE_BOX,isFixedRange,MinCT,MaxCT>
  : public NonsquareBox<Scalar_t,Eigen::Dynamic,DVO>
{
  private:
    static_assert(isFixedRange==false,
    "Compile-time box range for non-square box is not supported");
    using Base_t = NonsquareBox<Scalar_t,Eigen::Dynamic,DVO>;
  public:
    using Var_t = typename Base_t ::Var_t;

    inline void setMin(const Var_t & v)
    {
      Base_t::setMin(v);
    }

    inline void setMax(const Var_t & v)
    {
      Base_t::setMax(v);
    }

    inline int dimensions() const
    {
      assert(this->_minV.size()==this->_maxV.size());
      return this->_minV.size();
    }

};



template<typename Scalar_t,
  int Dim,DoubleVectorOption DVO,
  BoxShape BS,
  bool isFixedRange,
  TemplateVal_t<Scalar_t> MinCT,
  TemplateVal_t<Scalar_t> MaxCT>
using BoxDims = typename std::conditional<Dim==Eigen::Dynamic,
  BoxRTDim<Scalar_t,DVO,BS,isFixedRange,MinCT,MaxCT>,
  BoxCTDim<Scalar_t,Dim,DVO,BS,isFixedRange,MinCT,MaxCT>>::type;

} //internal

} //namespace Eigen

#endif // EIGEN_HEU_BOXRTCTDIMS_HPP
