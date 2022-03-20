// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.



#ifndef Heu_TYPES_HPP
#define Heu_TYPES_HPP

#include <Eigen/Core>
#include <type_traits>
#include <array>
#include <vector>
#include "Enumerations.hpp"
#include "Constants.hpp"

namespace Heu {

template<typename scalar_t,int Dim>
using stdContainer =
    typename std::conditional<Dim==Eigen::Dynamic,
        std::vector<scalar_t>,
        std::array<scalar_t,Dim>>::type;

template<int Size>
using stdVecD_t = stdContainer<double,Size>;


template<typename scalar_t,int Size,DoubleVectorOption DVO>
using Container = typename std::conditional<
    (DVO!=DoubleVectorOption::Eigen),
    stdContainer<scalar_t,Size>,
    Eigen::Array<double,Size,1>>::type;

//template<DoubleVectorOption dvo,size_t Dim>
//using FitnessVec_t= Container<double,Dim,dvo>;

namespace internal 
{

template<int _ObjNum>
struct heu_initializeSize
{
    template<typename A>
    inline static void resize(A * v,size_t sz) {
        assert(v->size()==sz);
    }
};

template<>
struct heu_initializeSize<Eigen::Dynamic>
{
    template<typename A>
    inline static void resize(A * v,size_t sz) {
        v->resize(sz);
    }
};

}   //  internal

}   //  Heu

#endif // Heu_TYPES_HPP
