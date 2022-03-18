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

#include <type_traits>
#include <array>
#include <vector>
#include "Enumerations.hpp"
#include "Constants.hpp"

namespace Heu {

template<typename scalar_t,size_t Dim>
using stdContainer =
    typename std::conditional<Dim==Runtime,
        std::vector<scalar_t>,
        std::array<scalar_t,Dim>>::type;

template<typename scalar_t,size_t Dim>
using EigenContainer =
    typename std::conditional<Dim==Runtime,
        Eigen::Array<scalar_t,Eigen::Dynamic,1>,
        Eigen::Array<scalar_t,Dim,1>>::type;


template<size_t Size>
using stdVecD_t = stdContainer<double,Size>;


template<typename scalar_t,size_t Size,DoubleVectorOption DVO>
using Container = typename std::conditional<
    (DVO!=DoubleVectorOption::Eigen),
    stdContainer<scalar_t,Size>,
    EigenContainer<scalar_t,Size>>::type;

///Array type when using Eigen array(s)
template<size_t Size>
using EigenVecD_t = EigenContainer<double,Size>;

//template<DoubleVectorOption dvo,size_t Dim>
//using FitnessVec_t= Container<double,Dim,dvo>;

template<size_t _ObjNum>
struct initializeSize
{
    template<typename A>
    inline static void resize(A * v,size_t sz) {
        assert(v->size()==sz);
    }
};

template<>
struct initializeSize<Runtime>
{
    template<typename A>
    inline static void resize(A * v,size_t sz) {
        v->resize(sz);
    }
};

};

#endif // Heu_TYPES_HPP
